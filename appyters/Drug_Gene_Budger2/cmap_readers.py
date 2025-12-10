import pandas as pd 
import numpy as np
import hashlib
import polars as pl

def prepare_novartis_data(gene, URL = 'https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/novartis_de'):
    '''
    gene: gene symbol to retrieve
    URL: Novartis data storage location

    output: results dataframe from Novartis data
    '''
    try:
        novartis_de = pd.read_feather(f'{URL}/{gene}.f').set_index('index')
    except:
        # print(f'{gene} not found in Novartis')
        return None
     # format p-values
    novartis_de['log10adj.P.Val'] = novartis_de['P.Adj'].replace(0,1e-323).map(np.log10)*-1
    # rename logFC column for concordance with Ginkgo columns
    novartis_de.rename(columns={'LogFC':'logFC', 'P.Adj':'adj.P.Val'}, inplace=True)
    return novartis_de

def prepare_lincs_data(gene, URL='https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/lincs_de'):
    '''
    gene: gene symbol to retrieve
    URL: LINCS data storage location

    output: results dataframe from LINCS data
    '''
    try:
        lincs_de = pd.read_feather(f'{URL}/{gene}.f')
    except:
        # print(f'{gene} not found in LINCS')
        return None
     # format p-values
    lincs_de['log10adj.P.Val'] = lincs_de['adj.P.Val'].replace(0,1e-323).map(np.log10)*-1
    # remove CRISPR KO perturbations
    lincs_ko_perturbs = pd.read_csv('https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/lincs_ko_perturbs.txt', sep='\t')
    lincs_de = lincs_de[~lincs_de['Drug'].isin(lincs_ko_perturbs.cmap_name.to_list())]
    return lincs_de

def hash_bucket(gene, num_buckets=512):
    '''
    gene: Gene symbol
    num_buckets: number of hash buckets to create

    output: integer hash for gene name (between 0-n_buckets)
    '''
    return int(hashlib.md5(gene.encode()).hexdigest(),16) % num_buckets

def prepare_tahoe_data(df, gene):
    '''
    df: DataFrame retrieved from Tahoe gene bucket file
    gene: gene to filter dataframe

    output: results dataframe from Tahoe data filtered to gene
    '''
    tahoe_de = df[df['gene_name']==gene]
    if tahoe_de.shape[0] == 0:
        # print(f'{gene} not found in Tahoe')
        return None
    tahoe_de['log10adj.P.Val'] = tahoe_de['padj'].replace(0,1e-323).map(np.log10)*-1
    tahoe_de.rename(columns = {'log2FoldChange':'logFC', 'drug':'Drug', 'padj':'adj.P.Val', 'group':'Perturbation', 'gene_name':'Gene'}, inplace=True)
    tahoe_de['GeneDir'] = np.where(tahoe_de['UpReg']==1,'Up','Dn')
    return tahoe_de

def retrieve_tahoe_data(gene_set, URL='https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/tahoe_de'):
    '''
    gene_set: list of gene symbols

    output: dictionary of tahoe data for each gene in gene_set
    '''
    hash_dict = {}
    for g in gene_set:
        hash_dict[g] = str(hash_bucket(g))
    hash_dict_rev = {}
    for gene, hash in hash_dict.items():
        if hash in hash_dict_rev:
            hash_dict_rev[hash].append(gene)
        else:
            hash_dict_rev[hash]=[gene]
    tahoe_results = {}
    for hash, genes in hash_dict_rev.items():
        df = pd.read_parquet(f'{URL}/gene_bucket_{hash}.parquet', use_pandas_metadata=False)
        for g in genes:
            tahoe_data = prepare_tahoe_data(df, g)
            tahoe_results[g] = tahoe_data
    return tahoe_results

def prepare_ginkgo_data_dict(gene, cell_types, URL='https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/ginkgo_de'):
    '''
    gene: gene symbol to retrieve
    cell types: Ginkgo cell types in GDPx1 and GDPx2
    URL: Ginkgo data storage location

    output: dictionary with dataframes for all Ginkgo cell types
    '''
    try:
        df = pd.read_feather(f'{URL}/{gene}.f')
    except:
        # print('Gene not found in Ginkgo')
        return None
    cell_type_results = {}
    for k in cell_types:
        subset = df[df['Perturbation'].str.contains(k)]
        subset['log10adj.P.Val'] = subset['adj.P.Val'].replace(0,1e-323).map(np.log10)*-1
        subset = subset.drop('index', axis=1)
        cell_type_results[k] = subset
        
    return cell_type_results

def prepare_ginkgo_data_df(gene, cell_types, URL='https://appyters.maayanlab.cloud/storage/DrugRegulators_Appyter/ginkgo_de'):
    '''
    gene: gene symbol to retrieve
    cell types: Ginkgo cell types in GDPx1 and GDPx2
    URL: Ginkgo data storage location

    output: results dataframe for all Ginkgo cell types
    '''
    try:

        df = pd.read_feather(f'{URL}/{gene}.f')
    except:
        # print('Gene not found in Ginkgo')
        return None
    cell_type_results = {}
    for k in cell_types:
        subset = df[df['Perturbation'].str.contains(k)]
        subset['log10adj.P.Val'] = subset['adj.P.Val'].replace(0,1e-323).map(np.log10)*-1
        subset = subset.drop('index', axis=1)
        cell_type_results[k] = subset
    
    all_df = pd.concat(cell_type_results.values(), ignore_index=True)
        
    return all_df