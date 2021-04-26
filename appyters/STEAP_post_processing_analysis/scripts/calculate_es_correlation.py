# to allow importing to work correctly (in a dirty way)
import os
import sys
import inspect
from utils import get_esmu_data
filepath = os.path.abspath(inspect.getfile(inspect.currentframe()))
currentdir = os.path.dirname(filepath)
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import pandas as pd
import itertools
from scipy.stats import spearmanr
from pathlib import Path


def calculate_spearmanr(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Uses spearman correlation to calculate the ES gene correlation.
    """
    corr_list = []
    for x, y in itertools.combinations(dataframe.columns, 2):
        corr_frame = dataframe.loc[:, [x, y]].fillna(0).copy()
        # ES value should be > 0 in both celltypes
        corr_frame = corr_frame[(corr_frame > 0).all(1)]
        corr, pval = spearmanr(
            corr_frame.iloc[:, 0].values,
            corr_frame.iloc[:, 1].values
        )
        corr_list.append([x, y, corr, pval])
    df = pd.DataFrame(
        corr_list,
        columns=['celltypex', 'celltypey', 'corr', 'pval']
    )
    return df


def correct_pval_correlation(corr_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Corrects the pvalues from the ES gene correlation using Bonferroni.
    '''
    correct_df = corr_df.copy()
    n_test = correct_df.shape[0]
    correct_df['pval_bonferroni'] = correct_df['pval'] * n_test
    correct_df.loc[correct_df['pval_bonferroni'] > 1, 'pval_bonferroni'] = 1
    return correct_df


def calculate_es_corr(datasets: list[str]) -> pd.DataFrame:
    """
    Calculates the expression specificity (ES) gene correlation between
    all celltypes in the input datasets.

    Parameters
    ----------
    datasets : list[str]
        List of the names of datasets. These names shouls correspong to
        the .csv file in the esmu directory.

    Returns
    -------
    es_corr_df : pd.DataFrame
        Pandas dataframe containing the ES gene correlation and pvalue
        between the gwas phenotypes.
    """
    df_list = []
    scrna_dict = get_esmu_data()
    for dataset in datasets:
        df_esmu = scrna_dict[dataset]
        df_esmu.columns = [f'{dataset}, {ct}' for ct in df_esmu.columns]
        df_list.append(df_esmu)

    merged_es_df = pd.concat(df_list, join='outer', axis=1)
    merged_es_df.sort_index(axis=1, inplace=True)
    es_corr_df = calculate_spearmanr(merged_es_df.fillna(0))
    es_corr_df = correct_pval_correlation(es_corr_df)
    return es_corr_df


if __name__ == "__main__":
    df_all = pd.read_hdf('data/CELLECT_output/data.h5', 'df_all')
    datasets = df_all['specificity_id'].unique().tolist()
    es_corr_df = calculate_es_corr(datasets)
    es_corr_df.to_hdf('data/CELLECT_output/data.h5', key='es_corr_df')
