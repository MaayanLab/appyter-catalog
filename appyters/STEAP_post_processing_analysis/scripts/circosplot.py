# to allow importing to work correctly (in a dirty way)
import os
import sys
import inspect
filepath = os.path.abspath(inspect.getfile(inspect.currentframe()))
currentdir = os.path.dirname(filepath)
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import constants
from scripts.pyCircos import Gcircle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import re
import pandas as pd
from typing import Union


def get_links_and_nodes(
    corr_df: pd.DataFrame,
    pivot_corr_df: pd.DataFrame,
    pattern: str,
    sign_threshold: int
) -> tuple[pd.DataFrame, list[str]]:
    """
    Generates the links and nodes for the circosplot.
    """
    # calculate mean corr across methods
    links = pivot_corr_df['corr'].mean(axis=1).reset_index().rename(
        columns={0: 'corr'}
    ).copy()
    # only get corr if in gwas_group_dict
    links = links[links[['gwasx', 'gwasy']].apply(
        lambda x: x.str.contains(pattern)
    ).all(axis=1)]
    nodes = pd.concat(
        [links['gwasx'], links['gwasy']]
    ).drop_duplicates().to_list()
    pthres = (0.05/(corr_df.shape[0]/len(constants.METHODS)))  # bonferroni
    sign_index = pivot_corr_df['pval'][
        pivot_corr_df['pval'] < pthres
    ].dropna(thresh=sign_threshold).index
    links.set_index(['gwasx', 'gwasy'], inplace=True)
    sign_links = links[links.index.isin(sign_index)].reset_index()
    return sign_links, nodes


def get_gwas_links(
    pivot_corr_df: pd.DataFrame,
    pattern: str,
    gwas_name: str
) -> pd.DataFrame:
    """
    Finds all correlations between the input gwas_name and gwas part of
    the pattern.
    """
    gwas_links = pivot_corr_df[(pivot_corr_df['gwasx'] == gwas_name)
                               |
                               (pivot_corr_df['gwasy'] == gwas_name)].copy()
    gwas_links = gwas_links[(gwas_links['gwasx'].str.contains(pattern))
                            |
                            (gwas_links['gwasy'].str.contains(pattern))]
    gwas_links['gwasx'], gwas_links['gwasy'] = np.where(
        gwas_links['gwasy'] == gwas_name,
        (gwas_links['gwasx'], gwas_links['gwasy']),
        (gwas_links['gwasy'], gwas_links['gwasx'])
    )
    gwas_links.sort_values('gwasx', inplace=True)
    return gwas_links


def chord_plot_pandas(
    gcircle: Gcircle,
    row: pd.Series,
    gwas_group_dict: dict[str, list[str]],
    corr_limit: tuple[float, float],
    bottom: float,
):
    """
    Pandas funtion to create chord plot.
    """
    if row['corr'] >= corr_limit[1] or row['corr'] <= corr_limit[0]:
        nodes = []
        for i in ['x', 'y']:
            node = row[f"gwas{i}"]
            for k, v in gwas_group_dict.items():
                if bool(re.search('|'.join(v), node)):
                    node = f"{k}, {node}"
            nodes.append(node)
        gcircle.chord_plot(
            [nodes[0], 0, 0, bottom],
            [nodes[1], 0, 0, bottom],
            color=cmap(norm(row['corr'])),
            alpha=.3
        )


def annotation_layer(
    gcircle: Gcircle,
    links: pd.DataFrame,
    nodes: list[str],
    color_bar: str,
    bottom: float,
    axis_color: str = 'k',
):
    """
    Creates the outer anotation barplots.
    """
    step = 0.001
    theta = np.linspace(0, (2-(2*step))*np.pi, len(nodes))
    scale = 500
    original_bottom = bottom
    for m in constants.METHODS:
        corr_values = links[f'corr_{m}'] * scale  # scales beta so its visible

        gcircle.ax.bar(
            theta,
            corr_values,
            color=cmap(norm(corr_values)),
            edgecolor='none',
            alpha=.5,
            width=step*15,
            bottom=bottom+scale,
            zorder=9
        )
        x = np.linspace(0, 2*np.pi, len(nodes))
        y = np.array([bottom+scale]*len(x))
        gcircle.ax.plot(x, y, alpha=1, color=axis_color, lw=.5, zorder=8)
        for i_scale in [-1, 1, -0.5, 0.5]:
            gcircle.ax.plot(
                x,
                y+i_scale*scale,
                alpha=.75,
                color=axis_color,
                lw=.25,
                zorder=8
            )
        bottom += 2*scale

    color_sign = ['white', 'white', color_bar]
    links_m = links.copy()
    links_m['color'] = 'white'  # 'None'
    for i, c in enumerate(color_sign, 1):
        pval_links = links_m.loc[:, links_m.columns.str.contains('pval')]
        significant_links = pval_links.sum(1) == i
        links_m.loc[significant_links, 'color'] = c

    gcircle.ax.bar(
        theta,
        y-original_bottom+scale,
        color=links_m['color'],
        edgecolor='none',
        alpha=.5,
        width=step*np.pi*5,
        bottom=original_bottom,
        zorder=0
    )


def preprocessing(
    corr_df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    gwas_name: Union[str, None],
    sign_threshold: int,
) -> tuple[pd.DataFrame, list[str], Union[pd.DataFrame, None]]:
    """
    Return links, nodes and gwas_links. The links (correlations) and
    nodes are used for the inner chord plot. The gwas links are used
    for the annotation layers.
    """
    pattern = "|".join(
        [v for values in gwas_group_dict.values() for v in values]
    )
    pivot_corr_df = corr_df.pivot(index=['gwasx', 'gwasy'], columns='method')
    links, nodes = get_links_and_nodes(
        corr_df,
        pivot_corr_df,
        pattern,
        sign_threshold
    )
    gwas_links = None
    if gwas_name is not None:
        pthres = (0.05/(corr_df.shape[0]/len(constants.METHODS)))  # bonferroni
        # check significance of correlations
        pivot_corr_df.reset_index(inplace=True)
        pivot_corr_df['pval'] = pivot_corr_df['pval'] < pthres
        pivot_corr_df.columns = [' '.join(col).strip().replace(' ', '_')
                                 for col in pivot_corr_df.columns.values]
        gwas_links = get_gwas_links(pivot_corr_df, pattern, gwas_name)
    return links, nodes, gwas_links


def plot(
    corr_df: pd.DataFrame,
    gwas_group_dict: dict[str, list[str]],
    gwas_name: Union[str, None] = None,
    sign_threshold: int = len(constants.METHODS),
    corr_limit: tuple[float, float] = (0, 0),
    bottom: float = 1200,
    ylim: float = 5000,
    color_map: str = 'tab10',
    figsize: tuple[float, float] = (10, 10),
    color_bar: str = 'tab:green',
    save: bool = False,
    filename: str = 'circosplot.png'
):
    """
    Returns a pandas dataframe showing which cell types are shared
    in which gwas phenotype group.

    Parameters
    ----------
    corr_df : pd.DataFrame
        Pandas dataframe containing the correlation and pvalue between
        the gwas phenotypes.
    gwas_group_dict : dict[str, list[str]]
        Dictionary with the phenotype group name as key and (regex) keywords
        of the phenotypes in a list as values.
        The values will be used as nodes.
    gwas_name : str or None
        The name of the gwas used in the input dataframe which is to be plotted
        in the outer circles as bar plots (annotation layer).
        If set to None the annotation layers won't be generated.
    sign_threshold : int
        Only use celltype correlation significant in 'sign_threshold' methods
        (by default number of methods).
    corr_limit : typle[float,float]
        Set a limit on which correlation to plot in the inner chord plot.
        For example (-0.5,0.5) will only plot correlations <-.5 and >0.5.
        Nonetheless, only significant correlations will be shown.
    bottom : float
        The position (radius in polar coordinates) where the nodes will
        be plotted.
    ylim : float
        Set the ylim of the plot (radius in polar coordinates).
    color_map : str
        The colormap used for the nodes. The nodes will be colorcoded based
        on which key the belong to in the gwas_group_dict.
    figsize : tuple[float, float]
        Width, height in inches.
    color_bar : str
        The color of the highlighted areas in the annotation layer.
    save : bool
        Whether to save the circoplot (default False).
    filename : str
        Filename (and path) of the saved circosplot.
    """
    gwas_group_dict = {k:v for k,v in gwas_group_dict.items() if gwas_name not in v}
    global cmap
    global norm
    links, nodes, gwas_links = preprocessing(
        corr_df,
        gwas_group_dict,
        gwas_name,
        sign_threshold
    )
    plt.style.use('default')
    cmap = cm.get_cmap('bwr')
    norm = plt.Normalize(-1, 1)
    cmap_group = {group_name: cm.get_cmap(color_map)(i)
                  for i, group_name in enumerate(gwas_group_dict.keys())}

    # make nodes
    gcircle = Gcircle()
    for node in nodes:
        for k, v in gwas_group_dict.items():
            if re.search("|".join(v), node):
                gcircle.add_locus(
                    f"{k}, {node}",
                    2,
                    bottom=bottom,
                    linewidth=1,
                    interspace=0,
                    facecolor=cmap_group[k],
                    edgecolor=cmap_group[k]
                )

    gcircle.set_locus(figsize=figsize)  # Create figure object
    gcircle.ax.set_ylim(0, ylim)
    # make chords inside the circle
    links.apply(
        lambda row: chord_plot_pandas(
            gcircle,
            row,
            gwas_group_dict,
            corr_limit,
            bottom
        ),
        axis=1
    )
    # make the barplots around the circle
    if gwas_name is not None:
        annotation_layer(
            gcircle,
            gwas_links,
            nodes,
            color_bar,
            bottom,
            axis_color='k'
        )
    # add legend
    patches = [mpatches.Patch(color=v, label=k) for k, v in cmap_group.items()]
    gcircle.ax.legend(
        handles=patches,
        bbox_to_anchor=(1.2, 1),
        loc=1,
        frameon=False
    )
    if save:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
