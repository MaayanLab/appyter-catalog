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
import re
import upsetplot
import matplotlib.pyplot as plt
from typing import Union


def create_group(gwas: str, gwas_group_dict: dict[str, list[str]]) -> str:
    """
    Finds the input string in the gwas_group_dict values and returns
    the key if found.
    """
    reverse_dict = {
        v: k for k, v_list in gwas_group_dict.items() for v in v_list
    }
    for k, v in reverse_dict.items():
        if bool(re.search(k, gwas)):
            return v


def add_count_to_group_dict(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]]
) -> dict[str, list[str]]:
    """
    Adds the number of gwas phenotypes in the group to the
    group name (key in dictionary).
    """
    unqiue_gwas = df['gwas'].unique().tolist()
    reverse_dict = {
        v: k for k, v_list in gwas_group_dict.items() for v in v_list
    }
    gwas_count = {k: 0 for k in reverse_dict.keys()}
    for gwas in unqiue_gwas:
        for k, v in reverse_dict.items():
            if bool(re.search(k, gwas)):
                gwas_count[k] += 1

    gwas_group_dict_copy = {}
    for k, v_list in gwas_group_dict.items():
        tot_count = (sum([gwas_count[v] for v in v_list]))
        gwas_group_dict_copy[f'{k} (N={tot_count})'] = v_list
    return gwas_group_dict_copy


def make_df_sliced(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    sign_threshold: int
) -> pd.DataFrame:
    """
    Creates a dataframe which is used as input for the upsetplot
    and excel file.
    """
    pattern = '|'.join(
        [v for v_list in gwas_group_dict.values() for v in v_list]
    )
    df_sliced = df[(df['gwas'].str.contains(pattern))].copy()
    df_sliced = df_sliced[
        df_sliced[f'pvalue_{constants.PVAL_CORRECTION}'] <= 0.05
    ]
    df_sliced = df_sliced.groupby(
        ['gwas', 'specificity_id', 'annotation']
    ).size().reset_index().rename(columns={0: 'N_sign'})
    # significant in sign_threshold methods
    df_sliced = df_sliced[df_sliced['N_sign'] >= sign_threshold]
    df_sliced['group'] = df_sliced.apply(
        lambda x: create_group(x['gwas'], gwas_group_dict),
        axis=1
    )
    return df_sliced


def make_df_upset(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    sign_threshold: int,
    sort_categories_by: Union[str, None]
) -> pd.DataFrame:
    """
    Creates a dataframe which is used as input for the upsetplot.
    """
    df_sliced = make_df_sliced(df, gwas_group_dict, sign_threshold)
    df_sliced_upset = df_sliced.groupby(
        ['specificity_id', 'annotation']
    )['group'].agg(list).reset_index()
    df_sliced_upset['group'] = df_sliced_upset['group'].apply(set)
    df_sliced_upset = df_sliced_upset['group'].value_counts().reset_index()
#     groups = list(df_sliced['group'].unique())
    groups = list(gwas_group_dict.keys())[::-1]
    for g in groups:
        df_sliced_upset[g] = df_sliced_upset['index'].apply(lambda x: g in x)

    df_sliced_upset.drop(columns='index', inplace=True)
    if sort_categories_by is None:
        new_column = groups.copy()
        new_column.append('group')
        df_sliced_upset = df_sliced_upset.reindex(columns=new_column)

    df_sliced_upset.set_index(groups, inplace=True)
    return df_sliced_upset


def get_shared_celltypes(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    sign_threshold: int = len(constants.METHODS)-1,
    save_to_excel: bool = False,
    filename: str = 'upsetplot.xlsx'
) -> pd.DataFrame:
    """
    Returns a pandas dataframe showing whihc cell types are shared in
    which gwas phenotype group.

    Parameters
    ----------
    df : pd.DataFrame
        The input pandas dataframe contains information about the phenotype,
        celltypes, method and enrichment (beta) values with corresponding
        p-values.
    gwas_group_dict : dict[str, list[str]]
        Dictionary with the phenotype group name as key and
        (regex) keywords
        of the phenotypes in a list as values.
    sign_threshold : int
        Only use celltypes significant in 'sign_threshold' methods
        (by default number of methods - 1).
    save_to_excel : bool
        Whether to save to an excel file (default False).
    filename : str
        Filename (and path) of the saved upsetplot.

    Returns
    -------
    dataframe : pd.DataFrame
        Output pandas dataframe showing the celltypes shared in the input
        gwas_group_dict gwas phenotype groups.
    """
    df_sliced = make_df_sliced(df, gwas_group_dict, sign_threshold)
    df_sliced = df_sliced.groupby(
        ['specificity_id', 'annotation']
    )[['group', 'gwas']].agg(set).reset_index()
    df_sliced['gwas'] = df_sliced['gwas'].apply(lambda x: ', '.join(list(x)))
    df_sliced['N_groups'] = df_sliced['group'].apply(len)
    df_sliced['group'] = df_sliced['group'].apply(lambda x: ', '.join(list(x)))
    df_sliced.sort_values(
        ['N_groups', 'group'],
        ascending=[False, True],
        inplace=True
    )
    df_sliced.drop(columns='N_groups', inplace=True)
    if save_to_excel:
        df_sliced.to_excel(filename, index=False)
    return df_sliced


def plot_upset(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    sign_threshold: int = len(constants.METHODS)-1,
    save: bool = False,
    filename: str = 'upsetplot.png',
    show_percentages: bool = False,
    sort_by: str = 'cardinality',
    with_lines: bool = True,
    element_size: Union[float, None] = 46,
    show_counts: Union[str, None] = '%d',
    sort_categories_by: Union[str, None] = 'cardinality'
):
    """
    Plots an upsetplot showing the shared number of celltypes between
    the different gwas phenotype groups.

    Parameters
    ----------
    df : pd.DataFrame
        The input pandas dataframe contains information about the
        phenotype, celltypes, method and enrichment (beta) values
        with corresponding p-values.
    gwas_group_dict : dict[str, list[str]]
        Dictionary with the phenotype group name as key and
        (regex) keywords of the phenotypes in a list as values.
    sign_threshold : int
        Only use celltypes significant in 'sign_threshold' methods
        (by default number of methods - 1).
    save : bool
        Whether to save plot (default False).
    filename : str
        Filename (and path) of the saved upsetplot.
    show_percentages : bool,
        Whether to label the intersection size bars with the percentage
        of the intersection relative to the total dataset.
        This may be applied with or without show_counts.
    sort_by : {'cardinality', 'degree'}
        If 'cardinality', subset are listed from largest to smallest.
        If 'degree', they are listed in order of the number of categories
        intersected.
    with_lines : bool
        Whether to show lines joining dots in the matrix, to mark multiple
        categories being intersected.
    element_size : float or None
        Side length in pt. If None, size is estimated to fit figure
    show_counts : bool or str,
        Whether to label the intersection size bars with the cardinality
        of the intersection. When a string, this formats the number.
        For example, '%d' is equivalent to True.
    sort_categories_by : {'cardinality', None}
        Whether to sort the categories by total cardinality, or leave them
        in the provided order.

    """
    gwas_group_dict_copy = add_count_to_group_dict(df, gwas_group_dict)
    df_sliced_upset = make_df_upset(
        df,
        gwas_group_dict_copy,
        sign_threshold,
        sort_categories_by
    )
    plt.style.use('default')
    upsetplot.plot(
        df_sliced_upset['group'],
        show_percentages=show_percentages,
        sort_by=sort_by,
        with_lines=with_lines,
        element_size=element_size,
        show_counts=show_counts,
        sort_categories_by=sort_categories_by
    )
    if save:
        plt.savefig(filename, dpi=200, bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    df_all = pd.read_hdf('data/data.h5', 'df_all')
    plot_upset(df_all, constants.GWAS_GROUP_DICT, save=True,
               filename='figures_and_tables/upsetplot.png')
    get_shared_celltypes(df_all, constants.GWAS_GROUP_DICT, save_to_excel=True,
                         filename='figures_and_tables/upsetplot.xslx')
