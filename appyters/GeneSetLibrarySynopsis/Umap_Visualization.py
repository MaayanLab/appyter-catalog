
from bokeh.palettes import Category20
from bokeh.models import HoverTool, ColumnDataSource, RangeSlider
from bokeh.plotting import figure, show, save, output_file
from bokeh.io import output_notebook
from sklearn.feature_extraction.text import TfidfVectorizer
from maayanlab_bioinformatics.enrichment import enrich_crisp
import matplotlib as mpl
import pandas as pd
import numpy as np
import urllib.request
import ssl
from collections import OrderedDict
import ipywidgets
import scanpy as sc
import os
import anndata
ssl._create_default_https_context = ssl._create_unverified_context


#### UMAP Visualization Class########

# attributes:
# query_set- list of genes that are to be queries against, detected for overlap ||| Type: List of Strings
# gene_libraries- list of enrichr libraries whose gene sets are to be analyzed for overlap with query gene set ||| List of Strings
# sig_value significance for fisher exact test ||| Type: float
# self.dataset : pandas dataframe where each row


class NoResults(Exception):
    pass


class APIFailure(Exception):
    pass


class NotValidFile(Exception):
    pass


class UMAP_Visualization:

    def __init__(self, query_set=[], gene_libraries=[], sig_value=.05, gmt_files=[]):
        self.query_set = [gene.strip() for gene in query_set]
        self.gene_libraries = gene_libraries
        self.significant_value = sig_value
        self.term_library_map = {}
        self.dataset = OrderedDict()
        self.dataset.update(self.process_gmt_files(gmt_files))
        self.dataset.update(self.get_libraries(gene_libraries))

    def process_gmt_files(self, gmt_files):
        # gmt_file: a list of gmt files to read
        lib_dict = OrderedDict()
        for file in gmt_files:
            path = os.getcwd()+f'\static\{file}'
            with open(path) as f:
                lines = f.readlines()
            for line in lines:
                parsed_line = line.split('\t')
                term, library, genes = parsed_line[0], parsed_line[1], parsed_line[2:]
                self.term_library_map[term] = library
                # genes = [x.split(',')[0].strip() for x in genes.split('\t')]
                genes[-1] = genes[-1][:-2]  # trim off newline characters!
                lib_dict[term] = ' '.join(genes)

        print(lib_dict)
        return lib_dict

    def get_libraries(self, libnames):
        lib_dict = OrderedDict()
        for libname in libnames:
            lib = self.get_scatter_libraries(libname)
            lib_dict.update(lib)

        return lib_dict

    def get_scatter_libraries(self, libname):
        print(f"\tOpening {libname} from Enrichr database...")
        with urllib.request.urlopen(f'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={libname}') as f:
            lines = [l.decode('utf-8') for l in f.readlines()]
        lib_dict = OrderedDict()
        for line in lines:
            tokens = line.split("\t\t")
            term = tokens[0]
            genes = [x.split(',')[0].strip() for x in tokens[1].split('\t')]
            lib_dict[term] = ' '.join(genes)
            # print(lib_dict[term])
            self.term_library_map[term] = libname
        print(lib_dict)
        return lib_dict

    def process_scatterplot(self, nneighbors=30, mindist=0.1, spread=1.0, maxdf=1.0, mindf=1):
        libdict = self.dataset
        print("\tTF-IDF vectorizing gene set data...")
        # computes tdfidf score--look this up
        vec = TfidfVectorizer(max_df=maxdf, min_df=mindf)
        X = vec.fit_transform(libdict.values())
        print(X.shape)
        adata = anndata.AnnData(X)
        adata.obs.index = libdict.keys()

        print("\tPerforming Leiden clustering...")
        # the n_neighbors and min_dist parameters can be altered
        sc.pp.neighbors(adata, n_neighbors=nneighbors)
        sc.tl.leiden(adata, resolution=1.0)
        sc.tl.umap(adata, min_dist=mindist, spread=spread, random_state=42)

        new_order = adata.obs.sort_values(by='leiden').index.tolist()
        adata = adata[new_order, :]
        adata.obs['leiden'] = 'Cluster ' + adata.obs['leiden'].astype('object')

        df = pd.DataFrame(adata.obsm['X_umap'])
        df.columns = ['x', 'y']

        df['cluster'] = adata.obs['leiden'].values
        df['term'] = adata.obs.index
        df['genes'] = [libdict[l] for l in df['term']]
        df['library'] = [self.term_library_map[l] for l in df['term']]

        return df

    def get_scatter_colors(self, df):
        clusters = pd.unique(df['cluster']).tolist()
        colors = list(Category20[20])[::2] + list(Category20[20])[1::2]
        color_mapper = {clusters[i]: colors[i % 20]
                        for i in range(len(clusters))}
        return color_mapper

    # def get_marker_mapper(self, df):
    #     markers = ["circle", "square", "triangle",
    #                "hex", "inverted_triangle", "diamond"]
    #     libs = pd.unique(df['library']).tolist()
    #     marker_mapper = {libs[i]: markers[i] for i in range(len(libs))}
    #     return marker_mapper

    def get_scatterplot(self, scatterdf):
        df = scatterdf.copy()
        color_mapper = self.get_scatter_colors(df)
        # marker_mapper = self.get_marker_mapper(df)
        df['color'] = df['cluster'].apply(lambda x: color_mapper[x])
        # df['marker'] = df['library'].apply(lambda x: marker_mapper[x])

        # range_slider = RangeSlider("title = Adjust x-axis",
        #                            start=0,
        #                            end=10,
        #                            step=1)

        tooltips = [
            ("Gene Set", "@gene_set"),
            ("Cluster", "@label"),
            ("Library", "@library")
        ]

        hover_emb = HoverTool(tooltips=tooltips)
        tools_emb = [hover_emb, 'pan', 'wheel_zoom', 'reset', 'save']

        plot_emb = figure(
            width=900,
            height=700,
            tools=tools_emb
        )

        source = ColumnDataSource(
            data=dict(
                x=df['x'],
                y=df['y'],
                gene_set=df['term'],
                colors=df['color'],
                label=df['cluster'],
                library=df['library'],
                # markers=df['marker']

            )
        )

        # hide axis labels and grid lines
        plot_emb.xaxis.major_tick_line_color = None
        plot_emb.xaxis.minor_tick_line_color = None
        plot_emb.yaxis.major_tick_line_color = None
        plot_emb.yaxis.minor_tick_line_color = None
        plot_emb.xaxis.major_label_text_font_size = '0pt'
        plot_emb.yaxis.major_label_text_font_size = '0pt'

        plot_emb.output_backend = "svg"

        plot_emb.xaxis.axis_label = "UMAP_1"
        plot_emb.yaxis.axis_label = "UMAP_2"

        s = plot_emb.scatter(
            'x',
            'y',
            size=4,
            source=source,
            color='colors',
            legend_group='label',
            # marker='markers'
        )

        plot_emb.add_layout(plot_emb.legend[0], 'right')

        return plot_emb
