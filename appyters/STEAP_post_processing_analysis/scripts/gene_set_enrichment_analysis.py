# to allow importing to work correctly (in a dirty way)
import os
import sys
import inspect
from utils import get_esmu_data, display_link
filepath = os.path.abspath(inspect.getfile(inspect.currentframe()))
currentdir = os.path.dirname(filepath)
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import constants
import pandas as pd
import mygene
import gseapy
from pathlib import Path
from typing import Union

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def get_annot_list(
    df: pd.DataFrame,
    name: str,
    param_list: list[str],
    rank: Union[int, None],
) -> list[list[str]]:
    """
    Returns a list of celltypes to be analyzed.
    """
    param = '|'.join(param_list)
    param_df = df[(df['gwas'].str.contains(param))]
    sign_enrichment = (param_df[f'pvalue_{constants.PVAL_CORRECTION}'] <= 0.05)
    sign_enrichment = sign_enrichment.astype(int)
    df_gsea = pd.concat(
        [
            param_df[['gwas', 'specificity_id', 'annotation']],
            sign_enrichment
        ], axis=1).rename(
        columns={f'pvalue_{constants.PVAL_CORRECTION}': 'count'}
    ).groupby(['gwas', 'specificity_id', 'annotation']).sum()
    n_gwas = df_gsea.index.unique(level='gwas').shape[0]
    # only use count==1 if multiple gwases are used in analysis
    if n_gwas == 1:
        # if only one gwas then include all celltypes
        min_count = 0
    else:
        # if multiple gwas then only cell types enriched in at least 2
        min_count = 1
#         min_count = 0

    df_gsea = df_gsea[df_gsea['count'] >= 2]  # only get if sign in >= methods
    df_gsea['count'] = 1  # reset count to 1
    df_gsea = df_gsea.groupby(
        ['specificity_id', 'annotation']
    ).sum().reset_index()
    N_total = df[(df['gwas'].str.contains(param))]['gwas'].unique().shape[0]
    df_gsea['freq'] = df_gsea['count'] / N_total
    df_gsea.sort_values('freq', ascending=False, inplace=True)
    df_gsea['rank'] = df_gsea['count'].rank(
        ascending=False,
        method='dense'
    ).astype(int)
    if rank is None:
        df_gsea = df_gsea[(df_gsea['count'] > min_count)]
    else:
        df_gsea = df_gsea[(df_gsea['rank'] <= rank)
                          &
                          (df_gsea['count'] > min_count)]

    annot_list = df_gsea[['specificity_id', 'annotation']].values.tolist()
    print(
        f"Top {rank} ranked cell-types (N={len(annot_list)}) in {name} GWAS:"
    )
    for row in df_gsea.iterrows():
        a = row[1]
#         print(f"{a[4]}. {a[0]}: {a[1]} (N={a[2]}, freq={a[3]:.2})")
        print(f"{a[4]}. {a[0]}: {a[1]} (N={a[2]})")

    print()
    return annot_list


def get_top_genes(annot_list: list[list[str]]) -> dict[str, list[str]]:
    """
    Returns the top genes specific to the celltype and converts
    the gene Enseble ID to Gene Symbol.
    """
    celltype_genes_dict = {}
    scrna_dict = get_esmu_data()
    for dataset, celltype in annot_list:
        df_esmu = scrna_dict[dataset]
        df = df_esmu[celltype]
        df = df.sort_values(ascending=False).head(
            round(constants.TOP_FREQ*len(df))
        )
        # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
        if df.shape[0] > 500 or df.shape[0] < 15:
            print(
                f'\033[93m{df.shape[0]} genes in {dataset}, {celltype} \
                which is outside the recommended range (15<G<500)\033[0m'
            )
        celltype_genes_dict[f'{dataset}, {celltype}'] = df.index.to_list()

    print('\nConverting Ensembl ID to Gene Symbol...')
    for i, (celltype, genes) in enumerate(celltype_genes_dict.items(), 1):
#         out_file = f"gsea_{celltype.replace(', ','-')}.xlsx"
#         if Path(out_file).is_file() and not constants.OVERWRITE_GSEA_ANALYSIS:
#             print(
#                 f'\n({i}/{len(celltype_genes_dict)}) \
#                 Converting {celltype} ({len(genes)} genes)'
#             )
#             print('Cell-type already analyzed. Skipping conversion...')
#         else:
#             print(
#                 f'\n({i}/{len(celltype_genes_dict)}) \
#                 Converting {celltype} ({len(genes)} genes)'
#             )
        with HiddenPrints():
            mg = mygene.MyGeneInfo()
            ginfo = mg.querymany(genes, scopes='ensembl.gene')
        gene_list = [g['symbol'] for g in ginfo if 'symbol' in g]
        celltype_genes_dict[celltype] = gene_list
    return celltype_genes_dict


def summarize_gsea(
    gsea_dict: dict[str, pd.DataFrame],
    correct_pval: bool = True,
    min_count: int = 0,
    save_to_excel: bool = False,
    filename: str = 'summarize_gsea.xlsx'
) -> pd.DataFrame:
    """
    Summarizes the gsea analysis by counting which term occurs
    the most often in the enriched cellypes.

    Parameters
    ----------
    gsea_dict : dict[str, pd.DataFrame]
        Dictionary where the keys are the analysed cell types and the values
        are a pandas dataframe containing significant gene-sets (terms)
        associated to the genes specific for the celltype.
    correct_pval : bool
        Corrects the pvalues for the number of celltypes analysed using
        the Bonferroni method.
    min_count : int
        Exclude terms which occur less or equal than 'min_counts'.
    save_to_excel : bool
        Whether to save to an excel file (default False).
    filename : str
        Filename (and path) of the saved file.

    Returns
    -------
    dataframe : pd.DataFrame
        Pandas dataframe summarizing the terms of the input celltypes by
        counting their occurences.
    """
    total = len(gsea_dict)
    df_list = []
    for k, v in gsea_dict.items():
        df = v.copy()
        df['Celltype'] = k
        df_list.append(df)

    gsea_df = pd.concat(df_list)
    if correct_pval:
        # Bonferroni correction
        gsea_df = gsea_df[(gsea_df['Adjusted P-value'] <= 0.05/total)]

    gsea_grouped_df = gsea_df.groupby(
        ['Gene_set', 'Term']
    )['Celltype'].agg(list).reset_index()
    col_name = f'Celltype_count (total={total})'
    gsea_grouped_df[col_name] = gsea_grouped_df['Celltype'].apply(
        lambda x: len(x)
    )
    gsea_grouped_df.sort_values(col_name, ascending=False, inplace=True)
    gsea_grouped_df['Celltype'] = gsea_grouped_df['Celltype'].apply(
        lambda x: '; '.join(x)
    )
    gsea_grouped_df = gsea_grouped_df[gsea_grouped_df[col_name] > min_count]
    if save_to_excel:
        gsea_grouped_df.to_excel(filename, index=False)
    return gsea_grouped_df


def gsea(
    df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    gene_set_list: list[str] = constants.GENE_SET_LIST,
    rank: Union[int, None] = constants.TOP_ANNOT
) -> dict[str, pd.DataFrame]:
    '''
    Performs gene-set enrichment analysis (GSEA) on enriched cell types.

    Parameters
    ----------
    df : pd.DataFrame
        The input pandas dataframe contains information about
        the phenotype, celltypes, method and enrichment (beta) values
        with corresponding p-values.
    gwas_group_dict : dict[str, list[str]]
        Dictionary with the phenotype group name as key and
        (regex) keywords of the phenotypes in a list as values.
    rank : int or None
        The minimum rank a cell type should have to be analysed.
        For example rank=5 means that only cell types ranked in
        the top 5 most occuring enriched cell types will be analysed.
        If rank=None, then all cell types will be analysed.

    Returns
    -------
    gsea_dict : dict[str, pd.DataFrame]
        Dictionary where the keys are the analysed cell types and
        the values are a pandas dataframe containing significant
        gene-sets (terms) associated to the genes specific for the celltype.
    '''
    print("Performing GSEA...\n")
    overwrite = constants.OVERWRITE_GSEA_ANALYSIS
    for name, param_list in gwas_group_dict.items():
        annot_list = get_annot_list(df, name, param_list, rank)
        celltype_genes_dict = get_top_genes(annot_list)

        print('\nRunning Enrichr...')
        gsea_dict = {}
        gsea_dir = 'gsea'
        for i, (celltype, genes) in enumerate(celltype_genes_dict.items(), 1):
            out_file = f"gsea_{celltype.replace(', ','-')}.xlsx"
            print(
                f'\n({i}/{len(celltype_genes_dict)}) \
                Analyzing {celltype} ({len(genes)} genes)...'
            )
            if Path(out_file).is_file() and not overwrite:
                print('Cell-type already analyzed. Skipping analysis...')
                gsea_dict[celltype] = pd.read_excel(out_file)
            else:
                df_list = []
                for gene_set in gene_set_list:
                    try:
                        enr = gseapy.enrichr(
                            gene_list=genes,
                            gene_sets=gene_set,
                            outdir=None,
                            organism='Human'
                        )
                        df = enr.results[
                            enr.results['Adjusted P-value'] <= 0.05
                        ]
                        df_list.append(df)
                    except Exception as e:
                        print(e)
                        continue
                try:
                    gsea_dict[celltype] = pd.concat(df_list)
                    gsea_dict[celltype].to_excel(out_file, index=False)
                    print('Writing to file...')
                except ValueError:
                    continue
    return gsea_dict


if __name__ == "__main__":
    df_all = pd.read_hdf('data/data.h5', 'df_all')
    gsea(df_all, constants.GWAS_GROUP_DICT)
