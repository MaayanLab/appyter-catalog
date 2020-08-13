"""Class for reading, parsing, and downloading data from the Drugmonizome API.
   Adapted fromo harmonizome.py.
"""

import gzip
import json
import os
import logging

# Support for both Python2.X and 3.X.
# -----------------------------------------------------------------------------
try:
    import io
    from urllib.request import urlopen
    from urllib.error import HTTPError
    from urllib.parse import quote_plus
except ImportError:
    from StringIO import StringIO
    from urllib2 import urlopen, HTTPError
    from urllib import quote_plus

try:
    input_shim = raw_input
except NameError:
    # If `raw_input` throws a `NameError`, the user is using Python 2.X.
    input_shim = input

import pandas as pd
import numpy as np
from scipy.sparse import lil_matrix, isspmatrix
from itertools import takewhile, repeat
from functools import reduce

def parse_gmt(fn, row_sep='\n', col_sep='\t'):
    '''
    Parser for reading Drugmonizome data in gmt format (ragged tsv)
    Each row is a drug set corresponding to a term
    First column is term ID, second column is empty, remaining columns are drug IDs (as InChI keys)
    
    Returns:
        dict: maps term ID to sets of associated drugs
    '''
    terms_to_drugs = {}
    with open(fn, 'r', newline=row_sep) as fh:
        for line in fh:
            lh = line.strip().split('\t')
            assert lh[1] == '', 'unexpected input format'
            terms_to_drugs[lh[0]] = set(lh[2:])
    return terms_to_drugs

def parse_gmt_to_df(fn, row_sep='\n', col_sep='\t'):
    '''
    Parser for reading Drugmonizome data in gmt format (ragged tsv)
    Each row is a drug set corresponding to a term
    First column is term ID, second column is empty, remaining columns are drug IDs (as InChI keys)
    
    Returns:
        dataframe: indices are drug IDs, columns are term IDs, filled with 1 if association and 0 otherwise
    '''
    terms_to_drugs = parse_gmt(fn, row_sep=row_sep, col_sep=col_sep)
    all_drugs = reduce(lambda s1, s2: s1.union(s2), terms_to_drugs.values())
    df = pd.DataFrame(0, index=sorted(all_drugs), columns=sorted(terms_to_drugs))
    for term in terms_to_drugs:
        df.loc[terms_to_drugs[term], term] = 1
    return df

def parse_multiple(fns, row_sep='\n', col_sep='\t'):
    '''
    Reads multiple Drugmonizome datasets and joins dataframe by shared drug IDs
    Parser for reading Drugmonizome data in gmt format (ragged tsv)
    Each row is a drug set corresponding to a term
    First column is term ID, second column is empty, remaining columns are drug IDs (as InChI keys)
    
    Returns:
        dataframe: indices are drug IDs, columns are term IDs, filled with 1 if association and 0 otherwise
    '''
    # Load individual datasets
    df_attributes = [parse_gmt_to_df(fn) for fn in fns]
    
    # Assemble all attribute datasets
    if len(df_attributes) > 1:
        # Obtain merged dataframe with omics and target data
        df = reduce(
            lambda a, b: pd.merge( # Merge two dataframes item by item
                a, # left
                b, # right
                # Items with the same left and right index are merged
                left_index=True,
                right_index=True,
                how='outer', # Keep mis-matched index
        ),
        df_attributes,
    )
    else:
        df = df_attributes[0]
    return df

def get_matches_df(drugmonizome_metadata, hits):
    '''
    Matches a list of drug screen hits to the appropriate Drugmonizome metadata.
    Looks for a full match between the name of the hit and the Name or a
    synonym of a drug.
    
    Param:
     - drugmonizome_metadata: dataframe with Drugmonizome metadata
     - hits: list/set of str names of drug screen hits
    
    Returns:
     - dataframe containing metadata for drug screen hits
    '''
    # format the set of hits
    hits = set(hit.strip().lower() for hit in hits if len(hit.strip()) > 0)
    print('Number of hits queried: {}'.format(len(hits)))
    
    # make boolean array for drugs in Drugmonizome, where a match to a hit is True
    in_name = np.array([drug.lower() in hits for drug in drugmonizome_metadata['Name']])
    in_synonyms = np.array([any(drug.strip().lower() in hits for drug in synonyms)
                            if isinstance(synonyms, list) else False
                            for synonyms in drugmonizome_metadata['Synonyms']])
    
    print('Number of matches: {} / {}'.format(np.sum(np.logical_or(in_synonyms, in_name)), len(drugmonizome_metadata)))
    
    # find unmatched drugs
    hits_name_copy = set(hits)
    for drug in drugmonizome_metadata['Name']:
        if drug.lower() in hits_name_copy:
            hits_name_copy.remove(drug.lower())

    hits_syn_copy = set(hits)
    for synonyms in drugmonizome_metadata['Synonyms']:
        if isinstance(synonyms, str): 
            for drug in synonyms:
                if drug.strip().lower() in hits_syn_copy:
                    hits_syn_copy.remove(drug.strip().lower())

    missing_hits = hits_name_copy.intersection(hits_syn_copy)
    print('Missing in Drugmonizome ({}): {}'.format(len(missing_hits), missing_hits))
    
    # filter hits metadata
    dfhits = drugmonizome_metadata.loc[np.logical_or(in_synonyms, in_name)]
    print('Total shape: {}'.format(dfhits.shape))
          
    return dfhits

# Enumerables and constants
# -----------------------------------------------------------------------------

class Enum(set):
    """Simple Enum shim since Python 2.X does not have them.
    """

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


def json_from_url(url):
    """Returns API response after decoding and loading JSON.
    """
    response = urlopen(url)
    data = response.read().decode('utf-8')
    return json.loads(data)


VERSION = 'v1'
API_URL = 'http://maayanlab.cloud/drugmonizome/data-api/api'
METADATA_URL = 'https://maayanlab.cloud/drugmonizome/metadata-api/entities'

# # This config objects pulls the names of the datasets, their directories, and
# # the possible downloads from the API. This allows us to add new datasets and
# # downloads without breaking this file.
# config = json_from_url('https://maayanlab.cloud/drugmonizome/metadata-api/libraries')
# DATASET_TO_LINK = {x['meta']['Library_name']: x['meta']['Download_link'] for x in config}

# Static version of Drugmonizome datasets, as the Drugmonizome resource is temporarily unavailable
LINK_BASE = 'https://appyters.maayanlab.cloud/storage/Drugmonizome_ML/DrugmonizomeGMTs/'
DATASET_TO_LINK = {
    'L1000FWD Downregulated GO Biological Processes': LINK_BASE + 'L1000FWD/L1000FWD_GO_Biological_Processes_drugsetlibrary_down.gmt',
    'L1000FWD Downregulated GO Cellular Components': LINK_BASE + 'L1000FWD/L1000FWD_GO_Cellular_Component_drugsetlibrary_down.gmt',
    'L1000FWD Downregulated GO Molecular Function': LINK_BASE + 'L1000FWD/L1000FWD_GO_Molecular_Function_drugsetlibrary_down.gmt',
    'L1000FWD Downregulated KEGG Pathways': LINK_BASE + 'L1000FWD/L1000FWD_KEGG_Pathways_drugsetlibrary_down.gmt',
    'L1000FWD Downregulated Signatures': LINK_BASE + 'L1000FWD/L1000FWD_signature_drugsetlibrary_down.gmt',
    'L1000FWD Predicted Side Effects': LINK_BASE + 'L1000FWD/L1000FWD_predicted_side_effects.gmt',
    'L1000FWD Upregulated GO Biological Process': LINK_BASE + 'L1000FWD/L1000FWD_GO_Biological_Processes_drugsetlibrary_up.gmt',
    'L1000FWD Upregulated GO Cellular Components': LINK_BASE + 'L1000FWD/L1000FWD_GO_Cellular_Component_drugsetlibrary_up.gmt',
    'L1000FWD Upregulated GO Molecular Function': LINK_BASE + 'L1000FWD/L1000FWD_GO_Molecular_Function_drugsetlibrary_up.gmt',
    'L1000FWD Upregulated KEGG Pathways': LINK_BASE + 'L1000FWD/L1000FWD_KEGG_Pathways_drugsetlibrary_up.gmt',
    'L1000FWD Upregulated Signatures': LINK_BASE + 'L1000FWD/L1000FWD_signature_drugsetlibrary_up.gmt',
    'Downregulated CREEDS Signatures': LINK_BASE + 'CREEDS/CREEDS_signature_drugsetlibrary_down.gmt',
    'Upregulated CREEDS Signatures': LINK_BASE + 'CREEDS/CREEDS_signature_drugsetlibrary_up.gmt',
    'DrugCentral Targets': LINK_BASE + 'DrugCentral/DrugCentral_target_drugsetlibrary.gmt',
    'DrugRepurposingHub Drug Targets': LINK_BASE + 'DrugRepurposingHub/DrugRepurposingHub_target_drugsetlibrary.gmt',
    'Drugbank Small Molecule Carriers': LINK_BASE + 'Drugbank/Drugbank_smallmolecule_carrier_drugsetlibrary.gmt',
    'Drugbank Small Molecule Enzymes': LINK_BASE + 'Drugbank/Drugbank_smallmolecule_enzyme_drugsetlibrary.gmt',
    'Drugbank Small Molecule Targets': LINK_BASE + 'Drugbank/Drugbank_smallmolecule_target_drugsetlibrary.gmt',
    'Drugbank Small Molecule Transporters': LINK_BASE + 'Drugbank/Drugbank_smallmolecule_transporter_drugsetlibrary.gmt',
    'Geneshot Associated Genes': LINK_BASE + 'Geneshot/Geneshot_associated_drugsetlibrary.gmt',
    'Geneshot Predicted AutoRIF Genes': LINK_BASE + 'Geneshot/Geneshot_predicted_autorif_drugsetlibrary.gmt',
    'Geneshot Predicted Coexpression Genes': LINK_BASE + 'Geneshot/Geneshot_predicted_coexpression_drugsetlibrary.gmt',
    'Geneshot Predicted Enrichr Genes': LINK_BASE + 'Geneshot/Geneshot_predicted_enrichr_drugsetlibrary.gmt',
    'Geneshot Predicted GeneRIF Genes': LINK_BASE + 'Geneshot/Geneshot_predicted_generif_drugsetlibrary.gmt',
    'Geneshot Predicted Tagger Genes': LINK_BASE + 'Geneshot/Geneshot_predicted_tagger_drugsetlibrary.gmt',
    'KinomeScan Kinases': LINK_BASE + 'KinomeScan/KinomeScan_kinase_drugsetlibrary.gmt',
    'PharmGKB Single Nucleotide Polymorphisms': LINK_BASE + 'PharmGKB/PharmGKB_snp_drugsetlibrary.gmt',
    'STITCH Targets': LINK_BASE + 'STITCH/STITCH_genes_drugsetlibrary.gmt',
    'ATC Codes Drugsetlibrary': LINK_BASE + 'ATC/ATC_drugsetlibrary.gmt',
    'DrugRepurposingHub Mechanisms of Action': LINK_BASE + 'DrugRepurposingHub/DrugRepurposingHub_moa_drugsetlibrary.gmt',
    'PharmGKB OFFSIDES Side Effects': LINK_BASE + 'PharmGKB/PharmGKB_OFFSIDES_side_effects_drugsetlibrary.gmt',
    'SIDER Indications': LINK_BASE + 'SIDER/SIDER_indications_drugsetlibrary.gmt',
    'SIDER Side Effects': LINK_BASE + 'SIDER/SIDER_side_effects_drugsetlibrary.gmt',
    'RDKIT MACCS Chemical Fingerprints': LINK_BASE + 'RDKIT/RDKIT_maccs_fingerprints_drugsetlibrary_updated.gmt',
}

# Drugmonizome class
# -----------------------------------------------------------------------------

class Drugmonizome(object):

    __version__ = VERSION
    DATASETS = DATASET_TO_LINK.keys()
    drug_metadata = None

    @classmethod
    def get(cls, entity, name=None, start_at=None):
        """Returns a single entity or a list, depending on if a name is
        provided. If no name is provided and start_at is specified, returns a
        list starting at that cursor position.
        """
        if name:
            name = quote_plus(name)
            return _get_by_name(entity, name)
        if start_at is not None and type(start_at) is int:
            return _get_with_cursor(entity, start_at)
        url = '%s/%s/%s' % (API_URL, VERSION, entity)
        result = json_from_url(url)
        return result

    @classmethod
    def next(cls, response):
        """Returns the next set of entities based on a previous API response.
        """
        start_at = _get_next(response)
        entity = _get_entity(response)
        return cls.get(entity=entity, start_at=start_at)

    @classmethod
    def download(cls, datasets=None):
        """For each dataset, creates a directory and downloads files into it.
        """
        # Why not check `if not datasets`? Because in principle, a user could 
        # call `download([])`, which should download nothing, not everything.
        # Why might they do this? Imagine that the list of datasets is
        # dynamically generated in another user script.
        if datasets is None:
            datasets = cls.DATASETS
            warning = 'Warning: You are going to download all Harmonizome '\
                      'data. This is roughly 30GB. Do you accept?\n(Y/N) '
            resp = input_shim(warning)
            if resp.lower() != 'y':
                return

        for dataset in datasets:
            if dataset not in cls.DATASETS:
                msg = '"%s" is not a valid dataset name. Check the `DATASETS`'\
                      ' property for a complete list of names.' % dataset
                raise AttributeError(msg)
            if not os.path.exists(dataset):
                os.mkdir(dataset)
            
            url = DATASET_TO_LINK[dataset]

            try:
                response = urlopen(url)
            except HTTPError as e:
                # Not every dataset has all downloads.
                raise Exception('Error downloading from %s: %s' % (url, e))
            
            filename = '%s/%s' % (dataset, url.split('/')[-1])

            if response.code != 200:
                raise Exception('This should not happen')
            
            if os.path.isfile(filename):
                logging.info('Using cached `%s`' % (filename))
            else:
                _download_file(response, filename)

            yield filename

    @classmethod
    def download_df(cls, datasets=None, **kwargs):
        for file in cls.download(datasets):
            yield _read_as_dataframe(file, **kwargs)
    
    @classmethod
    def get_datasets(cls):
        return cls.DATASETS
    
    @classmethod
    def read_drug_metadata(cls):
        """Reads all drug metadata into a dataframe
        """
        if cls.drug_metadata is None:
            entities = json_from_url(METADATA_URL)
            rows_list = []
            for entity in entities:
                rows_list.append(entity['meta'])
            cls.drug_metadata = pd.DataFrame(rows_list)
        return cls.drug_metadata

    @classmethod
    def get_InChI_keys(cls, hits):
        """Given list of drug names, finds matching InChI keys in Drugmonizome
        """
        df_drugs = cls.read_drug_metadata()
        df_hits = get_matches_df(df_drugs, hits)
        return list(df_hits['InChI_key'])
    
    @classmethod
    def map_names_to_inchi_keys(cls, hits):
        """Given list of drug names, finds matching InChI keys in Drugmonizome
           and returns dictionary mapping names to associated InChI keys
        """
        hits = set(hit.lower() for hit in hits)
        df_drugs = cls.read_drug_metadata()
        name_to_inchis = {}
        for i in range(len(df_drugs)):
            compound = df_drugs.iloc[i]
            if isinstance(compound.Synonyms, list):
                names = [compound.Name] + compound.Synonyms
            else:
                names = [compound.Name]
            for name in names:
                name = name.lower().strip()
                if name not in name_to_inchis:
                    name_to_inchis[name] = set()
                name_to_inchis[name].add(compound.InChI_key)
        return {hit: name_to_inchis[hit] for hit in hits if hit in name_to_inchis}
    
    @classmethod
    def get_drug_names(cls, inchi_keys):
        """Given list of InChI keys, finds matching drug names in Drugmonizome
        """
        df_drugs = cls.read_drug_metadata()
        df_drugs = df_drugs.set_index('InChI_key')
        return list(df_drugs.reindex(inchi_keys)['Name'])

# Utility functions
# -------------------------------------------------------------------------

def _get_with_cursor(entity, start_at):
    """Returns a list of entities based on cursor position.
    """
    url = '%s/%s/%s?cursor=%s' % (API_URL, VERSION, entity,str(start_at))
    result = json_from_url(url)
    return result


def _get_by_name(entity, name):
    """Returns a single entity based on name.
    """
    url = '%s/%s/%s/%s' % (API_URL, VERSION, entity, name)
    return json_from_url(url)


def _get_entity(response):
    """Returns the entity from an API response.
    """
    path = response['next'].split('?')[0]
    return path.split('/')[3]


def _get_next(response):
    """Returns the next property from an API response.
    """
    if response['next']:
        return int(response['next'].split('=')[1])
    return None


# This function was adopted from here: http://stackoverflow.com/a/15353312.
# def _download_and_decompress_file(response, filename):
#     """Downloads and decompresses a single file from a response object.
#     """
#     compressed_file = StringIO()
#     compressed_file.write(response.read())
#     compressed_file.seek(0)
#     decompressed_file = gzip.GzipFile(fileobj=compressed_file, mode='rb')
#     with open(filename, 'w+') as outfile:
#         outfile.write(decompressed_file.read())

def _download_file(response, filename):
    """
    """
    file = io.BytesIO(response.read())

    with open(filename, 'wb+') as outfile:
        outfile.write(file.read())


def json_ind_no_slash(ind_names, ind):
  return (
      json.dumps([ind_name.replace('/', '|')
                 for ind_name in ind_names]),
      [json.dumps([ii.replace('/', '|')
                  for ii in i])
       for i in ind],
  )

def _read_as_dataframe(fn):
    ''' Standard loading of dataframe '''
    # return fn
    import pandas as pd
    if fn.endswith('.gmt'):
        return parse_gmt_to_df(fn)
    else:
        raise Exception('Unable to parse this file into a dataframe.')
