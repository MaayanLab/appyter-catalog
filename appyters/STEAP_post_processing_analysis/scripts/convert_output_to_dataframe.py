# to allow importing to work correctly (in a dirty way)
import os
import sys
import inspect
filepath = os.path.abspath(inspect.getfile(inspect.currentframe()))
currentdir = os.path.dirname(filepath)
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import constants
import pandas as pd
from pathlib import Path
from statsmodels.stats.multitest import multipletests


def find_csv_file(directory: str) -> dict[str, dict[str, str]]:
    """
    Finds the priorirization.csv file in the input directory and uses the
    directory names to determine which method was used. This information is
    saced in a dictionary as:
    {phenotype : {method : path/to/file}}
    """
    file_dict = {}
    for path in Path(directory).rglob('prioritization.csv'):
        full_path = str(path)
        (name, method, __, __) = path.parts[-4:]
        name = name[8:]  # 'CELLECT-' is 8 long
        method = method[8:]
        if name not in file_dict:
            file_dict[name] = {}
        file_dict[name].update({method: full_path})
    return file_dict


def make_df(directory: str) -> pd.DataFrame:
    """
    Converts the prioritization.csv files in the directories to a
    pandas dataframe.

    Parameters
    ----------
    directory : str
        The output directory of CELLECT.
        There should be three subdirectories inside this directory.
        These subdirecories start with the name "CELLECT-" and
        end with the method used (H-MAGMA, LDSC or MAGMA)

    Returns
    -------
    dataframe : pd.DataFrame
        The output pandas dataframe contains information about the
        phenotype, celltypes, method and enrichment (beta) values
        with corresponding p-values.
    """
    file_dict = find_csv_file(directory)
    df_list_1 = []
    for name, d in file_dict.items():
        df_list_2 = []
        for method, file_path in d.items():
            df = pd.read_csv(file_path)
            df['method'] = method
            df.sort_values(
                by=['gwas', 'specificity_id', 'annotation'],
                inplace=True)
            df_list_2.append(df)
        df_list_1.extend(df_list_2)

    df_all = pd.concat(df_list_1, ignore_index=True)
    # count the number of methods (not used atm)
    df_all = df_all.merge(
        df_all.groupby(
            ['gwas', 'specificity_id', 'annotation']
        ).size().to_frame('n_methods'),
        on=['gwas', 'specificity_id', 'annotation'], how='left')
    # count the number of annotations/celltypes
    df_all.sort_values(by=['gwas', 'method'], inplace=True)
    df_all.reset_index(inplace=True, drop=True)
    return df_all


def pvalue_correction(
    dataframe: pd.DataFrame,
    method: str = 'bonferroni'
) -> pd.DataFrame:
    '''
    Corrects the pvalues in the input pandas dataframe for the
    multiple testing problem. The resulting output dataframe is
    the input dataframe with an additional corrected pvalues column.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The input pandas dataframe contains information about the phenotype,
        celltypes, method and enrichment (beta) values with corresponding
        p-values.
    method : str
        The pvalue correction method (default bonferroni).
        Other available methods are documented in
        https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    -------
    dataframe : pd.DataFrame
        The output pandas dataframe equals the input dataframe with an
        additional column containing the corrected pvalues.
    '''
    df_p = dataframe.pivot_table(
        values='pvalue',
        index=['method', 'gwas', 'specificity_id'],
        columns=['annotation']
    )
    df_p = df_p.apply(
        lambda row: multipletests(row.dropna(), method=method)[1],
        axis=1,
        result_type='reduce'
    ).apply(pd.Series).stack().reset_index().drop('level_3', axis=1)
    df_p.rename(columns={0: f"pvalue_{method}"}, inplace=True)
    df_p['annotation'] = dataframe['annotation']
    corrected_df = pd.merge(
        dataframe,
        df_p,
        on=['gwas', 'specificity_id', 'annotation', 'method']
    )
    return corrected_df


if __name__ == "__main__":
    df_all = make_df(constants.CELLECT_OUTDIR)
    df_all = pvalue_correction(df_all, method=constants.PVAL_CORRECTION)
    df_all.to_hdf('data/data.h5', key='df_all', mode='w')
