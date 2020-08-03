from PubChemQuery import PubChemQuery
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import requests

class DrugNameConverter:
    """
    Class for converting drug names to InChI keys using PubChem API to query drug names and RDKit for generating InChI keys.
    Includes options for using isomeric forms and for removing salts from drugs.
    """
    remover = SaltRemover()

    @classmethod
    def to_inchi_keys(cls, name, isomeric=True, strip_salts=True):
        """
        Queries PubChem API for a drug with a given name and returns a set of corresponding InChI Keys.

        Parameters:
            name (str):         name of drug
        
        Keyword arguments:
            isomeric (bool):    if True, returns InChI Keys computed from isomeric SMILES
                                otherwise, returns InChI Keys computed from canonical SMILES
            strip_salts (bool): if True, computed InChI Keys using both the original drug SMILES
                                and also the SMILES where all salts are removed
        
        Returns:
            inchi_keys (set): set of InChI Keys corresponding to the drug name queried
        """
        inchi_keys = set()
        for smiles in PubChemQuery.name_to_smiles(name, isomeric=isomeric):
            mol = Chem.MolFromSmiles(smiles)
            inchi_keys.add(Chem.MolToInchiKey(mol))
            if strip_salts:
                stripped_mol = cls.remover.StripMol(mol, dontRemoveEverything=True)
                inchi_keys.add(Chem.MolToInchiKey(stripped_mol))
        return inchi_keys

    @classmethod
    def batch_to_inchi_keys_single_thread(cls, names, verbose=0, **kwargs):
        """
        Queries PubChem API for a list of drug names and returns a dictionary mapping each name
        to a set of corresponding InChI Keys.

        Parameters:
            names (list or set):    drug names to query
        
        Keyword arguments:
            verbose (bool):         print progess if True
            **kwargs:               keyword arguments passed to cls.to_inchi_keys
        
        Returns:
            all_inchi_keys (dict):  dictionary mapping each drug name to a set of corresponding InChI Keys
        """
        all_inchi_keys = {}
        names = set(names)
        for name in names:
            inchi_keys = cls.to_inchi_keys(name, **kwargs)
            all_inchi_keys[name] = inchi_keys

            if verbose:
                print(f'Completed { len(all_inchi_keys) }/{ len(names) } drugs...', end='\r')
        return all_inchi_keys

    @classmethod
    def batch_to_inchi_keys(cls, names, num_cores=3, verbose=1, **kwargs):
        """
        Queries PubChem API for a list of drug names and returns a dictionary mapping each name
        to a set of corresponding InChI Keys. Uses multi-threading to parallelize requests.

        Parameters:
            names (list or set):    drug names to query
        
        Keyword arguments:
            num_cores (int/None):   number of threads to use; if None, uses CPU count (at least 1 and at most 12)
            verbose (bool):         show status bar if True
            **kwargs:               keyword arguments passed to cls.to_inchi_keys
        
        Returns:
            all_inchi_keys (dict):  dictionary mapping each drug name to a set of corresponding InChI Keys
        """
        requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound')  # initial request necessary before pooling (gives status code 400)
        if num_cores is None:
            num_cores = min(max(mp.cpu_count(), 1), 12)     # uses at least 1 core and at most 12
        names = list(set(names))
        with Pool(num_cores) as p:
            if verbose:
                res = list(tqdm(p.imap(partial(cls.to_inchi_keys, **kwargs), names), total=len(names)))
            else:
                res = p.map(partial(cls.to_inchi_keys, **kwargs), names)
        return dict(zip(names, res))
    
    @staticmethod
    def invert_dict(key_to_value_set: dict) -> dict:
        """
        Converts a dictionary with keys mapping to sets of values (e.g. drug name to set of InChI keys)
        into a dictionary with the values as keys, mapping to sets of the former keys (e.g. InChI key to drug names).
        """
        assert(isinstance(key_to_value_set, dict))
        value_to_key_set = {}
        for key in key_to_value_set:
            assert(isinstance(key_to_value_set[key], set))
            for value in key_to_value_set[key]:
                if value not in value_to_key_set:
                    value_to_key_set[value] = {key}
                else:
                    value_to_key_set[value].add(key)
        return value_to_key_set
