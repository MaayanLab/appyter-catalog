import pandas as pd
import numpy as np
import random
from tqdm import tqdm
import warnings
import itertools
import os
import sys

# utils file
import cmap_readers

def fetch_cmap_data(input_set):
    '''
    input_set: input gene list (iterable)

    output: dictionary of gene-specific dataframes, concatenated across all resources
    '''
    ginkgo_cell_types = ['human_epithelial_melanocytes', 'human_dermal_fibroblast', 'human_aortic_smooth_muscle_cells', 'human_skeletal_muscle_myoblasts','A549']
    novartis_perturbs  = {g:[] for g in input_set}
    # lincs_perturbs  = {g:[] for g in input_set}
    ginkgo_perturbs = {g:[] for g in input_set}

    not_found = {}
    for g in tqdm(input_set):
        novartis_perturbs[g] = cmap_readers.prepare_novartis_data(g)
        # lincs_perturbs[g] = cmap_readers.prepare_lincs_data(g)
        ginkgo_perturbs[g] = cmap_readers.prepare_ginkgo_data_df(g, ginkgo_cell_types)

    tahoe_perturbs = cmap_readers.retrieve_tahoe_data(input_set)
    # create single df for all resources
    cols_to_keep = ['Gene','Perturbation','Drug','logFC','GeneDir','PctRank', 'log10adj.P.Val']
    dataset_names = ['tahoe','novartis','ginkgo'] #'lincs'
    all_genes = {gene:None for gene in input_set}
    for g in all_genes.keys():
        all_df = pd.DataFrame()
        for i,pdict in enumerate([tahoe_perturbs, novartis_perturbs, ginkgo_perturbs]): #lincs_perturbs
            if pdict[g] is None:
                # print(f'{g} not found in {dataset_names[i]}')
                if dataset_names[i] in not_found:
                    not_found[dataset_names[i]].append(g)
                else:
                    not_found[dataset_names[i]] = [g]
                continue
            df = pdict[g][cols_to_keep]
            df['Dataset'] = dataset_names[i]
            if i == 0:
                all_df = df
            else:
                all_df = pd.concat([all_df, df], ignore_index=True)
        all_genes[g] = all_df
    return all_genes, not_found

def clean_drugname(d):
    if d == None:
        return d
    else:
        return d.casefold()
    
def add_pubchem_col(df, cid_df):
    '''
    df: single-gene dataframe
    cid_df: pubchem CID dataframe

    output: dataframe with CID column
    '''
    df['DrugMatch'] = df['Drug'].apply(clean_drugname)
    df = df.merge(cid_df, on='DrugMatch', how='left')
    df['DrugGroup'] = df['CID'].fillna(df['DrugMatch'])
    df = df.drop(columns=['DrugMatch'])
    return df

def rank_drugs(df, direction, pubchem_ids):
    '''
    df: CMAP dataframe 
    direction: up or down regulation

    output: dataframe ranked by up/down regulating signatures
    '''
    tmp = df.copy()
    # filter by direction
    tmp = tmp[tmp['GeneDir']==direction]
    # add drug matching column for downstream multi-gene analysis
    tmp = add_pubchem_col(tmp, pubchem_ids)
    # rank drugs by adjusted p-value
    tmp = tmp.sort_values(['log10adj.P.Val','PctRank'], ascending=[False,True], ignore_index=True)
    tmp['IndexRank'] = tmp.index + 1
    tmp['Rank'] = tmp['IndexRank'].rank(pct=True, ascending=True)
    return tmp.sort_values('Rank', ascending=True)

def agg_drug_rankings(df):
    '''
    df: single gene dataframe with signature rankings for up or down regulation

    output: dataframe with signature rankings averaged for each drug
    '''
    agg_rankings = df.groupby('DrugGroup').agg(
        MeanRank=('Rank','mean'), 
        NSigs=('Perturbation','size'),
        NResources = ('Dataset', 'nunique')) \
            .sort_values('MeanRank').reset_index()
    
    drug_match = df[['Drug','DrugGroup']].drop_duplicates(subset='DrugGroup', keep='first', ignore_index=True)
    agg_rankings = agg_rankings.merge(drug_match, on='DrugGroup', how='left')
    return agg_rankings

def single_gene_rank(gene_cmap_dict, input_set, pubchem_ids):
    '''
    gene_cmap_dict: dictionary from fetch_cmap
    input_set: input gene set
    pubchem_ids: dataframe of pubchem IDs

    Rank drugs for all genes in input set individually

    out: dictionary of drug ranks for individual genes
    '''
    # prepare drug rank output object
    drug_ranks = {g:{dir:None for dir in ['Up','Dn']} for g in input_set}
    # genes
    remove = []
    for g in drug_ranks.keys():
        # no CMAP resources had data for the gene
        if gene_cmap_dict[g].empty:
            remove.append(g)
            continue
        # direction
        for dir in drug_ranks[g].keys():
            rank_df = rank_drugs(gene_cmap_dict[g], direction=dir, pubchem_ids=pubchem_ids)
            drug_ranks[g][dir]= agg_drug_rankings(rank_df)
    # remove genes for which there was no CMAP data
    drug_ranks = {g:data for g,data in drug_ranks.items() if g not in remove}
    return drug_ranks

def combine_gene_results (agg_rank_dict, direction):
    '''
    agg_rank_dict: dictionary with drug rankings for multiple genes
    direction: 'Up' or 'Dn' regulation

    output: dataframe with drug rankings for given direction and method for multiple genes
    '''
    results = {'gene':[],
               'DrugGroup':[],
               'Drug':[],
               'MeanRank':[]}
    for g in agg_rank_dict.keys():
        results['gene'].extend([g] * agg_rank_dict[g][direction].shape[0])
        results['DrugGroup'].extend(agg_rank_dict[g][direction]['DrugGroup'].to_list())
        results['Drug'].extend(agg_rank_dict[g][direction]['Drug'].to_list())
        results['MeanRank'].extend(agg_rank_dict[g][direction]['MeanRank'].to_list())
    return pd.DataFrame(results)

def get_multi_gene_ranks (agg_drug_ranks, pubchem_ids):
    '''
    combined_df: dataframe from combine_gene_results

    output: dataframe of ranked drugs, scores averaged across genes
    '''
    multi_gene_ranks = {direction: None for direction in ['Up','Dn']}
    combined_gene_results = {direction: None for direction in ['Up','Dn']}

    for dir in ['Up','Dn']:
        combined_df = combine_gene_results(agg_drug_ranks, dir)
        mg_ranks = combined_df.groupby(['DrugGroup']).agg(MultiGeneRank = ('MeanRank','mean'), NGenes = ('gene', 'nunique')).reset_index().sort_values('MultiGeneRank', ascending=True)
        mg_ranks = mg_ranks.merge(pubchem_ids, left_on='DrugGroup',right_on='CID',how='left').drop(columns='CID')
        combined_df = combined_df.merge(pubchem_ids, left_on='DrugGroup',right_on='CID',how='left').drop(columns='CID')
        multi_gene_ranks[dir] = mg_ranks
        combined_gene_results[dir] = combined_df
    
    return combined_gene_results, multi_gene_ranks