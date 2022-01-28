# Basic libraries
import pandas as pd
import requests, json
import time
import numpy as np
import warnings
import random

# Visualization
import seaborn as sns
import scipy.stats as ss
import IPython
from IPython.display import HTML, display, Markdown, IFrame, FileLink
from itertools import combinations
from scipy import stats


# Data analysis
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform
from sklearn import cluster
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
import umap
from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
import scanpy as sc
import anndata
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
from maayanlab_bioinformatics.dge.limma_voom import limma_voom_differential_expression
from maayanlab_bioinformatics.enrichment.crisp import enrich_crisp, fisher_overlap
from statsmodels.stats.multitest import multipletests
from scipy.stats.mstats import gmean
# Bokeh
from bokeh.io import output_notebook
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, CustomJS, ColumnDataSource, Span, Select, Legend, PreText, Paragraph, LinearColorMapper, ColorBar, CategoricalColorMapper
from bokeh.layouts import layout, row, column, gridplot
from bokeh.palettes import all_palettes
import colorcet as cc
from bokeh.palettes import Category20

sc.settings.verbosity = 0

def check_files(fname):
    if fname == "":
        raise IOError
    if fname.endswith(".txt") == False and fname.endswith(".csv") ==False and fname.endswith(".tsv")==False:
        raise IOError
def check_df(df, col):
    if col not in df.columns:
        raise IOError



def display_statistics(data, description=""):
    print(description)
    print("Sample size:", data.n_obs)
    print("Feature size:", data.n_vars)
    
def load_seurat_files(mtx_filename, gene_filename, barcodes_filename, bool_adt=False):
    # bool_adt: if data contains both RNA and ADT 
    
    adata = anndata.read_mtx(mtx_filename).T
    with open(barcodes_filename, "r") as f:
        cells = f.readlines()
        cells = [x.strip() for x in cells]
    genes = pd.read_csv(
        gene_filename,
        header=None,
        sep='\t',
    )

    adata.var['gene_ids'] = genes.iloc[:, 0].values  
    if genes.shape[1] > 1:
        adata.var['gene_symbols'] = genes.iloc[:, 1].values
    else:
        adata.var['gene_symbols'] = genes.iloc[:, 0].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    adata.obs['barcode'] = cells
    adata.obs_names = cells
    adata.obs_names_make_unique(join="-")
    
    
    if bool_adt == True:
        genes.columns = ["gene_ids", "gene_symbols", "gene_annot"]
        adt_list = genes.loc[genes["gene_annot"]=="Antibody Capture", "gene_symbols"].tolist()
        adata_adt = adata[:, adata.var["gene_symbols"].isin(adt_list)]
        adata = adata[:, ~adata.var["gene_symbols"].isin(adt_list)]
        return adata, adata_adt
    else:
        return adata

def load_metadata(adata, meta_data_filename, meta_class_column_name):
    if meta_data_filename is not None and meta_data_filename != "":
        if meta_data_filename.endswith(".csv"):
            meta_df = pd.read_csv(meta_data_filename, index_col=0)
        else:
            meta_df = pd.read_csv(meta_data_filename, sep="\t", index_col=0)
        if meta_class_column_name == "":
            raise Exception ("Run time error: Please provide a proper column name for sample classes in metadata")
        try:
            check_df(meta_df, meta_class_column_name)
        except:
            raise Exception (f"Error! Column '{meta_class_column_name}' is not in metadata")
        adata.obs[meta_class_column_name] = meta_df.loc[:, meta_class_column_name]
        adata.var_names_make_unique()

    return adata

def load_data(dataset_name, rnaseq_data_filename, adt_data_filename, mtx_data_filename_rna, gene_data_filename_rna, barcode_data_filename_rna, mtx_data_filename_adt, gene_data_filename_adt, barcode_data_filename_adt, meta_data_filename=None, meta_class_column_name=None, table_counter=1):
    adata = None
    adata_adt = None
    
    # plain text
    if rnaseq_data_filename != "":
        display(Markdown(f"### Loading...{dataset_name}"))
        check_files(rnaseq_data_filename)
        check_files(adt_data_filename)
        # load rna
        try:
            if rnaseq_data_filename.endswith(".csv"):
                expr_df = pd.read_csv(rnaseq_data_filename, index_col=0).sort_index()
            else:
                expr_df = pd.read_csv(rnaseq_data_filename, index_col=0, sep="\t").sort_index()

            # convert df into anndata
            # adata matrix: sample x gene
            adata = anndata.AnnData(expr_df.T)
            adata.X = adata.X.astype('float64')    

        except:
            print("Error! Input files are in a wrong format. \
            Please check if the index of the expression data are genes and the columns are sample IDs. \
            Sample IDs in the expression data and the metadata should be matched")
        del expr_df
        
        # load adt
        try:
            if adt_data_filename.endswith(".csv"):
                expr_df = pd.read_csv(adt_data_filename, index_col=0).sort_index()
            else:
                expr_df = pd.read_csv(adt_data_filename, index_col=0, sep="\t").sort_index()

            # convert df into anndata
            # adata matrix: sample x gene
            adata_adt = anndata.AnnData(expr_df.T)
            adata_adt.X = adata_adt.X.astype('float64')    

        except:
            print("Error! Input files are in a wrong format. \
            Please check if the index of the expression data are genes and the columns are sample IDs. \
            Sample IDs in the expression data and the metadata should be matched")

        del expr_df
    # mtx files
    elif mtx_data_filename_rna != "":
        display(Markdown(f"### Loading...{dataset_name}"))
        if mtx_data_filename_adt == "":
            adata, adata_adt = load_seurat_files(mtx_data_filename_rna, gene_data_filename_rna, barcode_data_filename_rna, bool_adt=True)
        else:
            adata = load_seurat_files(mtx_data_filename_rna, gene_data_filename_rna, barcode_data_filename_rna)
            adata_adt = load_seurat_files(mtx_data_filename_adt, gene_data_filename_adt, barcode_data_filename_adt)
    
    if adata is not None:
        # load metadata
        adata = load_metadata(adata, meta_data_filename, meta_class_column_name)
        adata_adt = load_metadata(adata_adt, meta_data_filename, meta_class_column_name)
        
        # add batch info
        if meta_class_column_name == "":
            meta_class_column_name = "batch"
        adata.obs["batch"] = dataset_name
        adata_adt.obs["batch"] = dataset_name
        adata.obs.index = adata.obs.index + "-" + dataset_name
        adata_adt.obs.index = adata_adt.obs.index + "-" + dataset_name
        
        # common samples
        common_samples = list(set(adata_adt.obs.index.tolist()).intersection(adata.obs.index.tolist()))
        if len(common_samples) == 0:
            raise Exception("There are no matched samples.")
        adata = adata[common_samples, :]
        adata_adt = adata_adt[common_samples, :]
        
        
        table_counter = display_object(table_counter, f"Raw RNA data of {dataset_name}. The table displays the first 5 rows of the quantified RNA-seq expression dataset. Rows represent genes, columns represent samples, and values show the number of mapped reads.", adata.to_df().iloc[:10,:5].T.head(), istable=True)
        table_counter = display_object(table_counter, f"Raw protein data of {dataset_name}. The table displays the first 5 rows of the protein dataset. Rows represent genes, columns represent samples, and values show the number of mapped reads.", adata_adt.to_df().iloc[:10,:5].T.head(), istable=True)
        table_counter = display_object(table_counter, f"Metadata in {dataset_name}. The table displays the metadata associated with the samples in the RNA-seq dataset. Rows represent RNA-seq samples, columns represent metadata categories.", adata.obs.head(), istable=True)
        table_counter = display_object(table_counter, f"Sample size for each class in {dataset_name}. The table displays the number of samples in each class.", adata.obs.reset_index().groupby(meta_class_column_name).count(), istable=True)
        display_statistics(adata, f"### Statistics of RNA data in {dataset_name} ###")
        display_statistics(adata_adt, f"### Statistics of protein data in {dataset_name} ###")
    return adata, adata_adt, table_counter, meta_class_column_name
    
    

def create_download_link(df, title = "Download CSV file: {}", filename = "data.csv"):  
    if filename.endswith(".csv"):
        df.to_csv(filename)
    elif filename.endswith(".h5ad"): #anndata
        df.write(filename)
    html = "<a href=\"./{}\" target='_blank'>{}</a>".format(filename, title.format(filename))
    return HTML(html)
    
def display_link(url, title=None):
    if title is None:
        title = url
    raw_html = '<a href="%s" target="_blank">%s</a>' % (url, title)
   
    return display(HTML(raw_html))

def display_object(counter, caption, df=None, istable=True):
    if df is not None:
        display(df)
    if istable == True:
        display(Markdown("*Table {}. {}*".format(counter, caption)))
    else:
        display(Markdown("*Figure {}. {}*".format(counter, caption)))
    counter += 1
    return counter


def autoselect_color_by(sample_metadata):
    '''Automatically select a column in the sample_metadata df for coloring.
    '''
    color_by = None
    color_type = 'categorical'
    meta_col_nuniques = sample_metadata.nunique()
    # pick a column with the cardinality between 2 and 10
    meta_col_nuniques = meta_col_nuniques.loc[meta_col_nuniques.between(1, 30)]
    if len(meta_col_nuniques) > 0:
        color_by = meta_col_nuniques.index[0]
    else: # pick a numeric column
        is_number = np.vectorize(lambda x: np.issubdtype(x, np.number))
        meta_col_dtypes = sample_metadata.dtypes
        try:
            meta_col_is_number = is_number(meta_col_dtypes)
            if meta_col_is_number.sum() > 0:
                color_by = meta_col_dtypes.loc[meta_col_is_number].index[0]
                color_type = 'continuous'
        except:
            pass       
        
    return color_by, color_type

import hashlib
def str_to_int(string, mod):
    byte_string = bytearray(string, "utf8")
    return int(hashlib.sha256(byte_string).hexdigest(), base=16)%mod

def clr(x):    
    return np.log(x) - np.log(gmean(x))

def normalize(adata, normalization_method, log_normalization):
    tmp_adata = adata.copy()
    if normalization_method == "Seurat":
        sc.pp.filter_cells(tmp_adata, min_genes=200)
        sc.pp.filter_genes(tmp_adata, min_cells=3)
        sc.pp.normalize_total(tmp_adata, target_sum=1e4)        
        if log_normalization:
            sc.pp.log1p(tmp_adata)
        sc.pp.scale(tmp_adata, max_value=10)
    elif normalization_method == "CLR":
        tmp_adata.X = clr((tmp_adata.to_df()+1).values)
    return tmp_adata
    
    
def run_clustergrammer(dataset, meta_class_column_name, magic_normalization=False, nr_genes=800, metadata_cols=None, filter_samples=True,gene_list=None):
    # Subset the expression DataFrame using top 800 genes with largest variance
    if magic_normalization == True:
        data = dataset.uns["magic"]
    else:
        data = dataset.to_df().T
    meta_df = dataset.obs
    variances = np.var(data, axis=1)
    srt_idx = variances.argsort()[::-1]
    if gene_list == None or len(gene_list) == 0:
        expr_df_sub = data.iloc[srt_idx].iloc[:nr_genes]
    else:
        gene_list = gene_list.split("\n")
        common_gene_list = list(set(gene_list).intersection(set(data.index)))
        expr_df_sub = data.loc[common_gene_list, :]
        assert len(expr_df_sub.index) > 0
    # prettify sample names
    sample_names = ['::'.join([y, x]) for x,y in
                       zip(meta_df[meta_class_column_name], expr_df_sub.columns)]
    expr_df_sub.columns = sample_names
    expr_df_sub.index = ["Gene: "+str(x) for x in expr_df_sub.index]
    sample_name = ["Sample: "+x for x in sample_names]
    expr_df_sub.columns = sample_name


    treatment_type = ["Class: "+ x.split("::")[1] for x in sample_names]
    new_series = pd.DataFrame(treatment_type).T
    new_series.columns = expr_df_sub.columns
    expr_df_sub = pd.concat([new_series, expr_df_sub], axis=0)

    index_list = list(expr_df_sub.index)
    index_list = ["" if "Gene" not in str(x) else x for x in index_list]
    expr_df_sub.index = index_list
    #subset of expr_df_sub
    if len(expr_df_sub.columns) > 50:
        print("Input data is too large. Random sampling (n=50) is performed.")
        expr_df_sub = expr_df_sub.sample(50, axis=1)
    expr_df_sub_file = "expr_df_sub_file.txt"
    expr_df_sub.to_csv("expr_df_sub_file.txt", sep='\t')
    # POST the expression matrix to Clustergrammer and get the URL
    clustergrammer_url = 'https://maayanlab.cloud/clustergrammer/matrix_upload/'
    r = requests.post(clustergrammer_url, files={'file': open(expr_df_sub_file, 'rb')}).text

    return r
    
#############################################
########## 2. Plot
#############################################

def plot_clustergrammar(clustergrammer_url):
    clustergrammer_url = clustergrammer_url.replace("http:", "https:")
    display_link(clustergrammer_url, clustergrammer_url)
    # Embed
    display(IPython.display.IFrame(clustergrammer_url, width="1000", height="1000"))

def plot_scatter(umap_df, values_dict, option_list, sample_names, caption_text, category_list_dict=None, location='right', category=True, dropdown=False, figure_counter=0, x_axis_label="UMAP_1", y_axis_label="UMAP_2"):
    
    # init plot
    source = ColumnDataSource(data=dict(x=umap_df["x"], y=umap_df["y"], values=values_dict[option_list[0]], names=sample_names))
    if location == 'right':
        plot = figure(plot_width=800, plot_height=600)   
    else:
        plot = figure(plot_width=600, plot_height=600+20*len(category_list_dict[option_list[0]]))   
    if category == True:
        unique_category_dict = dict()
        for option in option_list:
            unique_category_dict[option] = sorted(list(set(values_dict[option])))
        
        # map category to color
        # color is mapped by its category name 
        # if a color is used by other categories, use another color
        factors_dict = dict()
        colors_dict = dict()
        for key in values_dict.keys():
            unused_color = list(Category20[20])
            factors_dict[key] = category_list_dict[key]
            colors_dict[key] = list()
            for category_name in factors_dict[key]:
                color_for_category = Category20[20][str_to_int(category_name, 20)]
                
                if color_for_category not in unused_color:
                    if len(unused_color) > 0:
                        color_for_category = unused_color[0]                        
                    else:
                        color_for_category = Category20[20][19]
                
                colors_dict[key].append(color_for_category)
                if color_for_category in unused_color:
                    unused_color.remove(color_for_category)
                    
        color_mapper = CategoricalColorMapper(factors=factors_dict[option_list[0]], palette=colors_dict[option_list[0]])
        legend = Legend()
        
        plot.add_layout(legend, location)
        scatter = plot.scatter('x', 'y', source=source, color={'field': 'values', 'transform': color_mapper}, legend_field="values")
        plot.legend.label_width = 30
        plot.legend.click_policy='hide'
        plot.legend.spacing = 1
        if location == 'below':
            location = 'bottom_left'
        plot.legend.location = location
        plot.legend.label_text_font_size = '10pt'
    else:
        color_mapper = LinearColorMapper(palette=cc.CET_D1A, low=min(values_dict[option_list[0]]), high=max(values_dict[option_list[0]]))
        color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12)
        plot.add_layout(color_bar, 'right')
        plot.scatter('x', 'y', source=source, color={'field': 'values', 'transform': color_mapper})
    
    tooltips = [
        ("Sample", "@names"),
        ("Value", "@values"),
    ]
    plot.add_tools(HoverTool(tooltips=tooltips))
    plot.output_backend = "svg"
    
    plot.xaxis.axis_label = x_axis_label#"UMAP_1"
    plot.yaxis.axis_label = y_axis_label#"UMAP_2"
    plot.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    plot.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
    plot.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    plot.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
    plot.xaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    plot.yaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels

    default_text = "Figure {}. {}{}"
    pre = Paragraph(text = default_text.format(figure_counter, caption_text, option_list[0]), width=500, height=100, style={"font-family":'Helvetica', "font-style": "italic"})
    figure_counter += 1
    if dropdown == True:
        if category == True:
            callback_adt = CustomJS(args=dict(source=source, \
                                              pre=pre, \
                                              values_dict=values_dict, \
                                              figure_counter=figure_counter, \
                                              color_mapper=color_mapper,\
                                              unique_category_dict=unique_category_dict,\
                                              category_list_dict=category_list_dict,\
                                              factors_dict=factors_dict,\
                                              colors_dict=colors_dict,\
                                              plot=plot,\
                                              scatter=scatter,
                                              caption_text=caption_text
                                             ), code="""        
                const val = cb_obj.value;                    
                source.data.values = values_dict[val]    
                color_mapper.factors = category_list_dict[val]
                color_mapper.palette = colors_dict[val]
                plot.legend = unique_category_dict[val]
                pre.text = "Figure "+figure_counter+". "+caption_text+val+"."; 
                plot.height = 600+20*(category_list_dict[val].length)
                source.change.emit();
            """)
        else:
            callback_adt = CustomJS(args=dict(source=source, \
                                              pre=pre, \
                                              values_dict=values_dict, \
                                              figure_counter=figure_counter,
                                              caption_text=caption_text), code="""        
                const val = cb_obj.value;    
                source.data.values = values_dict[val]    
                pre.text = "Figure "+figure_counter+". "+caption_text+val+".";  
                source.change.emit();
            """)

        # init dropdown menu
        select = Select(title="Select an option:", value=option_list[0], options=option_list)
        select.js_on_change('value', callback_adt)

        col = column(select, row(column(plot, pre)))
        show(col)
    else:
        col = column(plot, pre)
        show(col)
    return figure_counter

def average_graphs(graphs, weights=None):
    """Take the maximum edge value from each graph.

    graphs : iterable of sparse matrices
    weights : iterable of floats with same length as graphs
    """
    if weights is None:
        weights = [1 for g in graphs]
    if len(graphs) != len(weights):
        raise ValueError('graphs and weights need to have the same lenghts')
    out = weights[0] * graphs[0].copy()
    for i in range(1, len(graphs)):
        out += weights[i] * graphs[i]
    out /= sum(weights)
    return out

def clustering(adata, adata_adt, dataset, bool_plot, figure_counter, batch_correction=True):
    # clustering
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.umap(adata, min_dist=0.1)
    
    sc.tl.leiden(adata_adt, resolution=1.0)
    sc.tl.umap(adata_adt, min_dist=0.1)
    
    # joint clustering   
    joint = adata.copy()
    joint.obsm['protein'] = adata_adt.to_df()

    joint.uns['neighbors'] = {}
    joint.uns['neighbors']['connectivities'] = average_graphs([adata.uns['neighbors']['connectivities'], adata_adt.uns['neighbors']['connectivities']])
    joint.uns['neighbors']['connectivities_key'] = 'connectivities'
    joint.uns['neighbors']['distances_key'] = 'distances'
    joint.uns['neighbors']['params'] = {'n_neighbors': 30, 'method': 'umap', 'random_state': 0, 'metric': 'euclidean'}
    
    if batch_correction== True:
        sc.external.pp.bbknn(joint, batch_key="batch")
    else:
        sc.pp.neighbors(joint, n_neighbors=30)
    sc.tl.leiden(joint, resolution=1.0)
    sc.tl.umap(joint, min_dist=0.1)
    
    
    # sort by clusters
    new_order = joint.obs.sort_values(by='leiden').index.tolist()
    adata = adata[new_order, :]
    adata_adt = adata_adt[new_order, :]
    joint = joint[new_order, :]
    
    if bool_plot == True:
        umap_df = pd.DataFrame(joint.obsm['X_umap'])
        umap_df.columns = ['x', 'y']

        values_dict = dict()
        values_dict["Cluster"] = joint.obs["leiden"].values
        category_list_dict = dict()
        category_list_dict["Cluster"] = list(sorted(joint.obs["leiden"].values.unique()))
        figure_counter = plot_scatter(umap_df, values_dict, ["Cluster"], joint.obs.index.tolist(), "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False, figure_counter=figure_counter)

        display(create_download_link(joint.obs["leiden"], filename=f"clustering_{dataset}.csv"))
    return joint, figure_counter
def protein_levels(adata, adata_adt, joint, bool_plot, figure_counter):
    new_order = joint.obs.sort_values(by='leiden').index.tolist()
    adata = adata[new_order, :]
    adata_adt = adata_adt[new_order, :]
    joint = joint[new_order, :]
    if bool_plot == True:
        umap_df = pd.DataFrame(joint.obsm['X_umap'])
        umap_df.columns = ['x', 'y']

        option_list = [x for x in adata_adt.to_df().columns.tolist()]    
        tmp_adata_adt_norm_selected = adata_adt[:, option_list]
        values_dict = dict(zip(tmp_adata_adt_norm_selected.to_df().T.index.tolist(), tmp_adata_adt_norm_selected.to_df().T.values))
        figure_counter = plot_scatter(umap_df, values_dict, option_list, adata_adt.obs.index.tolist(), "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category=False, dropdown=True, figure_counter=figure_counter)
    return figure_counter


def differential_gene_expression_analysis(joint, diff_gex_method, enrichment_groupby, table_counter):
    
    if enrichment_groupby == "user_defined_class":
        classes = joint.obs[meta_class_column_name].unique().tolist()
        bool_cluster=False
        if len(classes) < 2:
            print("Warning: Please provide at least 2 classes in the metadata")
    elif enrichment_groupby == "Cluster":
        meta_class_column_name = "leiden"
        classes = sorted(joint.obs["leiden"].unique().tolist())
        classes.sort(key=int)
        bool_cluster=True
    else:
        meta_class_column_name = enrichment_groupby
        classes = sorted(joint.obs[meta_class_column_name].unique().tolist())
        classes.sort()
        bool_cluster=True
        
    if len(classes) > 5 and joint.n_obs > 5000:
        if diff_gex_method == "wilcoxon":
            print('Warning: There are too many cells/clusters. It cannot execute the analysis code for the data. The appyter randomly select 5000 samples. If you want to execute it with the whole data, please run it locally.')
        else:
            print('Warning: There are too many cells/clusters. It cannot execute the analysis code for the data. The appyter switched to Wilcoxon rank-sum method and randomly select 5000 samples. If you want to execute it with the whole data, please run it locally.')
            diff_gex_method = "wilcoxon"
        # randomly select 5K samples
        random_selected_samples = random.sample(joint.obs.index.tolist(), 5000)
        joint_random_sampled = joint[random_selected_samples, :]
        signatures = get_signatures(classes, joint_random_sampled, method=diff_gex_method, meta_class_column_name=meta_class_column_name, cluster=bool_cluster)
    else:
        signatures = get_signatures(classes, joint, method=diff_gex_method, meta_class_column_name=meta_class_column_name, cluster=bool_cluster)
    if len(classes) > 1:
        for label, signature in signatures.items():
            table_counter = display_object(table_counter, f"Top 5 Differentially Expressed Genes in {label}", signature.head(5), istable=True)
            display(create_download_link(signature, filename="DEG_{}.csv".format(label)))
    return signatures, bool_cluster, table_counter

def visualize_enrichment_analysis(joint, signatures, meta_class_column_name, diff_gex_method, enrichr_libraries_filename, enrichr_libraries, enrichment_groupby, libraries_tab, gene_topk, bool_cluster, bool_plot, figure_counter, table_counter):
    
    if diff_gex_method == "characteristic_direction":
        fc_colname = "CD-coefficient"
        sort_genes_by = "CD-coefficient"
        ascending = False
    elif diff_gex_method == "limma":
        fc_colname = "logFC"
        sort_genes_by = "t"
        ascending = False
    elif diff_gex_method == "edgeR":
        fc_colname = "logFC"
        sort_genes_by = "PValue"
        ascending = True
    elif diff_gex_method == "DESeq2":
        fc_colname = "log2FoldChange"
        sort_genes_by = "padj"
        ascending = True
    elif diff_gex_method == "wilcoxon":
        fc_colname = "logfoldchanges"
        sort_genes_by = "pvals_adj"
        ascending = True
    results = {}    
    if libraries_tab == 'Yes' or libraries_tab == 'All':
        results['enrichr'] = {}
        for label, signature in signatures.items():
            # Run analysis
            if enrichment_groupby == "user_defined_class":
                case_name = label.split(" vs. ")[1]
                col_name = meta_class_column_name
            elif enrichment_groupby == "Cluster":
                case_name = label.split(" vs. ")[0]
                col_name = "leiden"
            else:
                case_name = label.split(" vs. ")[0]
                col_name = "batch"

            results['enrichr'][label] = run_enrichr(signature=signature, signature_label=label, fc_colname=fc_colname,geneset_size=gene_topk, sort_genes_by = sort_genes_by,ascending=ascending)
            display(Markdown("*Enrichment Analysis Result: {} (Up-regulated in {})*".format(label, case_name)))
            display_link("https://maayanlab.cloud/Enrichr/enrich?dataset={}".format(results['enrichr'][label]["upregulated"]["shortId"]))
            display(Markdown("*Enrichment Analysis Result: {} (Down-regulated in {})*".format(label, case_name)))
            display_link("https://maayanlab.cloud/Enrichr/enrich?dataset={}".format(results['enrichr'][label]["downregulated"]["shortId"]))
        table_counter = display_object(table_counter, "The table displays links to Enrichr containing the results of enrichment analyses generated by analyzing the up-regulated and down-regulated genes from a differential expression analysis. By clicking on these links, users can interactively explore and download the enrichment results from the Enrichr website.", istable=True)
    if libraries_tab == 'No' or libraries_tab == 'All':
        results['user_defined_enrichment'] = {}
        for label, signature in signatures.items():

            # Run analysis
            if enrichment_groupby == "user_defined_class":
                case_name = label.split(" vs. ")[1]
                col_name = meta_class_column_name
            elif enrichment_groupby == "Cluster":
                case_name = label.split(" vs. ")[0]
                col_name = "leiden"
            else:
                case_name = label.split(" vs. ")[1]
                col_name = "batch"

            user_library = library_to_dict(get_library(enrichr_libraries_filename))

            # Sort signature
            up_signature = signature[signature[fc_colname] > 0].sort_values(sort_genes_by, ascending=ascending)
            up_genes = [x.upper() for x in up_signature.index[:gene_topk].tolist()]

            results['user_defined_enrichment'][label] = dict()
            _, _, results['user_defined_enrichment'][label]['enrichment_dataframe'] = enrichment_analysis(up_genes, user_library)
            results['user_defined_enrichment'][label]['enrichment_dataframe']['gene_set_library'] = enrichr_libraries_filename
    
    if "Gene Ontology" in enrichr_libraries:
        results['go_enrichment'] = {}
        for label, signature in signatures.items():
            # Run analysis
            results['go_enrichment'][label] = get_enrichr_results_by_library(results['enrichr'][label], label, library_type='go', version='2018')
            
    if "Pathway" in enrichr_libraries:
        # Initialize results
        results['pathway_enrichment'] = {}

        # Loop through results
        for label, signature in signatures.items():
            # Run analysis
            results['pathway_enrichment'][label] = get_enrichr_results_by_library(results['enrichr'][label], label, library_type='pathway')
    if "Transcription Factor" in enrichr_libraries:
        # Initialize results
        results['tf_enrichment'] = {}
        # Loop through results
        for label, signature in signatures.items():
            # Run analysis
            results['tf_enrichment'][label] = get_enrichr_result_tables_by_library(enrichr_results=results['enrichr'][label], signature_label=label, library_type='tf')
    if "Kinase" in enrichr_libraries:
        # Initialize results
        results['kinase_enrichment'] = {}

        # Loop through results
        for label, enrichr_results in results['enrichr'].items():
            # Run analysis
            results['kinase_enrichment'][label] = get_enrichr_result_tables_by_library(enrichr_results=enrichr_results, signature_label=label, library_type="ke")

    if "miRNA" in enrichr_libraries:
        results['mirna_enrichment'] = {}

        # Loop through results
        for label, enrichr_results in results['enrichr'].items():
            # Run analysis
            results['mirna_enrichment'][label] = get_enrichr_result_tables_by_library(enrichr_results=enrichr_results, signature_label=label, library_type="mirna")
    if "Cell Type" in enrichr_libraries:
        results['celltype_enrichment'] = {}
        for label, signature in signatures.items():
            # Run analysis
            results['celltype_enrichment'][label] = get_enrichr_results_by_library(results['enrichr'][label], label, library_type='celltype')
    if "Disease" in enrichr_libraries:
        results['disease_enrichment'] = {}
        for label, signature in signatures.items():
            # Run analysis
            results['disease_enrichment'][label] = get_enrichr_results_by_library(results['enrichr'][label], label, library_type='disease', version='2018')
    
    library_option_list = set()
    for label, signature in signatures.items():
        if bool_cluster == True:
            cluster_names = [label.split(" vs. ")[0].replace("Cluster ", "")]
        else:
            cluster_names = label.split(" vs. ")

        for key in results.keys():
            if key.endswith("enrichment") == False:
                continue
            enrichment_results = results[key][label]
            meta_df = joint.obs
            if bool_cluster == True:
                for cluster_name in cluster_names:
                    for direction in ['upregulated']:
                        if direction in enrichment_results:
                            enrichment_dataframe = enrichment_results[direction]
                        else:
                            enrichment_dataframe = enrichment_results["enrichment_dataframe"]
                            
                        if enrichment_dataframe.empty == True:
                            raise Exception("Enrichment analysis returns empty results. Please check if your data contains proper gene names.")
                        libraries = enrichment_dataframe['gene_set_library'].unique() 
                        for library in libraries:
                            enrichment_dataframe_library = enrichment_dataframe[enrichment_dataframe['gene_set_library']==library]
                            top_term = enrichment_dataframe_library.iloc[0]['term_name']
                            if library not in meta_df.columns:
                                meta_df.insert(0, library, np.nan)
                            meta_df[library] = meta_df[library].astype("object")
                            if bool_cluster == True:
                                meta_df.loc[meta_df[col_name]==cluster_name, library] = top_term
                            else:
                                meta_df.loc[meta_df[meta_class_column_name]==cluster_name, library] = top_term
                            library_option_list.add(library)
            else: # bool_cluster == False
                for direction in ['upregulated']:
                    if direction in enrichment_results:
                        enrichment_dataframe = enrichment_results[direction]
                    else:
                        enrichment_dataframe = enrichment_results["enrichment_dataframe"]
                    if enrichment_dataframe.empty == True:
                        raise Exception("Enrichment analysis returns empty results. Please check if your data and library contains shared items.")

                            
                    libraries = enrichment_dataframe['gene_set_library'].unique()  
                    for library in libraries:
                        enrichment_dataframe_library = enrichment_dataframe[enrichment_dataframe['gene_set_library']==library]
                        top_term = enrichment_dataframe_library.iloc[0]['term_name']
                        if library not in meta_df.columns:
                            meta_df.insert(0, library, np.nan)
                        meta_df[library] = meta_df[library].astype("object")
                        if direction == "upregulated":
                            meta_df.loc[meta_df[meta_class_column_name]==cluster_names[0], library] = top_term
                        else:
                            meta_df.loc[meta_df[meta_class_column_name]==cluster_names[1], library] = top_term
                        library_option_list.add(library)

    library_option_list = list(library_option_list)
    # umap info into dataframe 
    umap_df = pd.DataFrame(joint.obsm['X_umap'])
    umap_df.columns = ['x', 'y']

    option_list = library_option_list  
    adata_norm_selected = joint.obs[option_list].fillna("NaN")
    
    values_dict = dict(zip(adata_norm_selected.T.index.tolist(), adata_norm_selected.T.values))
    category_list_dict = dict()
    for option in option_list:
        category_list_dict[option] = list(sorted(adata_norm_selected[option].unique()))
    if bool_plot == True:
        figure_counter = plot_scatter(umap_df, values_dict, option_list, joint.obs.index.tolist(), "Scatter plot of the samples. Each dot represents a sample and it is colored by enriched terms in library ", location='below', category_list_dict=category_list_dict, category=True, dropdown=True, figure_counter=figure_counter)
    return joint, option_list, figure_counter, table_counter




def summary(joint, option_list, table_counter):
    for col in option_list:
        counts = joint.obs[['barcode', col]].groupby(col).count()
        counts.columns = ['# of Samples']
        counts["Percentage (%)"] = counts['# of Samples']/counts['# of Samples'].sum() * 100
        counts = counts.sort_values("Percentage (%)", ascending=False)
        table_counter = display_object(table_counter, "The number of samples for each category in {}".format(col), counts, istable=True)
    return table_counter

def get_signatures(classes, dataset, method, meta_class_column_name, cluster=True, filter_genes=True):
           
    robjects.r('''edgeR <- function(rawcount_dataframe, g1, g2) {
        # Load packages
        suppressMessages(require(limma))
        suppressMessages(require(edgeR))
        colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
        rownames(colData) <- c(g1,g2)
        colnames(colData) <- c("group")
        colData$group = relevel(as.factor(colData$group), "Control")
        y <- DGEList(counts=rawcount_dataframe, group=colData$group)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        et <- exactTest(y)
        res <- topTags(et, n=Inf)
        # Return
        res <- as.data.frame(res)
        results <- list("edgeR_dataframe"= res, "rownames"=rownames(res))
        return (results)
    }
    ''')

    robjects.r('''deseq2 <- function(rawcount_dataframe, g1, g2) {
        # Load packages
        suppressMessages(require(DESeq2))
        colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
        rownames(colData) <- c(g1,g2)
        colnames(colData) <- c("group")
        colData$group = relevel(as.factor(colData$group), "Control")
        dds <- DESeqDataSetFromMatrix(countData = rawcount_dataframe, colData = colData, design=~(group))
        dds <- DESeq(dds)
        res <- results(dds)
        res[which(is.na(res$padj)),] <- 1
        res <- as.data.frame(res)
        results <- list("DESeq_dataframe"= res, "rownames"=rownames(res))
        return(results)
    }
    ''')
    
    expr_df = dataset.to_df().T
    raw_expr_df = dataset.raw.to_adata().to_df().T
    meta_df = dataset.obs
    
    signatures = dict()

    sc.tl.rank_genes_groups(dataset, meta_class_column_name, method='t-test', use_raw=True)
    if cluster == True:
        
        # cluster 0 vs rest        
        for cls1 in classes:
            signature_label = " vs. ".join(["Cluster {}".format(cls1), "rest"])
            cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls1].index.tolist() #case
            non_cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]!=cls1].index.tolist() #control
            sample_ids = non_cls1_sample_ids.copy()
            sample_ids.extend(cls1_sample_ids)
            tmp_raw_expr_df = raw_expr_df[sample_ids]
            if method == "limma":
                signature = limma_voom_differential_expression(tmp_raw_expr_df.loc[:, non_cls1_sample_ids], tmp_raw_expr_df.loc[:, cls1_sample_ids])
            elif method == "characteristic_direction":
                signature = characteristic_direction(expr_df.loc[:, non_cls1_sample_ids], expr_df.loc[:, cls1_sample_ids], calculate_sig=False)
            elif method == "edgeR":
                edgeR = robjects.r['edgeR']
                edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(non_cls1_sample_ids), pandas2ri.conversion.py2rpy(cls1_sample_ids)))

                signature = pd.DataFrame(edgeR_results[0])
                signature.index = edgeR_results[1]
                signature = signature.sort_values("PValue", ascending=True)
            elif method == "DESeq2":
                DESeq2 = robjects.r['deseq2']
                DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(non_cls1_sample_ids), pandas2ri.conversion.py2rpy(cls1_sample_ids)))

                signature = pd.DataFrame(DESeq2_results[0])
                signature.index = DESeq2_results[1]
                signature = signature.sort_values("padj", ascending=True)
            elif method == "wilcoxon":   
                
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls1).set_index('names').sort_values('pvals', ascending=True)
                dedf = dedf.loc[dataset.var.index, :].sort_values("scores", ascending=False)
                signature = dedf
                
            signatures[signature_label] = signature
    else:
        for cls1, cls2 in combinations(classes, 2):
            signature_label = " vs. ".join([cls1, cls2])
            cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls1].index.tolist() #control
            cls2_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls2].index.tolist() #case
            sample_ids = cls1_sample_ids.copy()
            sample_ids.extend(cls2_sample_ids)
            tmp_raw_expr_df = raw_expr_df[sample_ids]
            if method == "limma":
                signature = limma_voom_differential_expression(tmp_raw_expr_df.loc[:, cls1_sample_ids], tmp_raw_expr_df.loc[:, cls2_sample_ids])
                
            elif method == "characteristic_direction":
                signature = characteristic_direction(expr_df.loc[:, cls1_sample_ids], expr_df.loc[:, cls2_sample_ids], calculate_sig=False)
            elif method == "edgeR":
                edgeR = robjects.r['edgeR']
                edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))

                signature = pd.DataFrame(edgeR_results[0])
                signature.index = edgeR_results[1]
                signature = signature.sort_values("PValue", ascending=True)
            elif method == "DESeq2":
                DESeq2 = robjects.r['deseq2']
                DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))

                signature = pd.DataFrame(DESeq2_results[0])
                signature.index = DESeq2_results[1]
                signature = signature.sort_values("padj", ascending=True)
            elif method == "wilcoxon":   
                sc.tl.rank_genes_groups(dataset, meta_class_column_name, method='wilcoxon', use_raw=True)
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls2).set_index('names').sort_values('pvals', ascending=True)
                dedf = dedf.loc[dataset.var.index, :].sort_values("scores", ascending=False)
                signature = dedf
            signatures[signature_label] = signature
    return signatures






def submit_enrichr_geneset(geneset, label):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(geneset)
    payload = {
        'list': (None, genes_str),
        'description': (None, label)
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    time.sleep(0.5)
    data = json.loads(response.text)
    return data


def run_enrichr(signature, signature_label, geneset_size=500, fc_colname = 'logFC', sort_genes_by='t', ascending=True):

    # Sort signature
    up_signature = signature[signature[fc_colname] > 0].sort_values(sort_genes_by, ascending=ascending)
    down_signature = signature[signature[fc_colname] < 0].sort_values(sort_genes_by, ascending=ascending)
    
    # Get genesets
    genesets = {
        'upregulated': up_signature.index[:geneset_size],
        'downregulated': down_signature.index[-geneset_size:]
    }

    # Submit to Enrichr
    enrichr_ids = {geneset_label: submit_enrichr_geneset(geneset=geneset, label=signature_label+', '+geneset_label+', from Bulk RNA-seq Appyter') for geneset_label, geneset in genesets.items()}
    enrichr_ids['signature_label'] = signature_label
    return enrichr_ids

def get_enrichr_results(user_list_id, gene_set_libraries, overlappingGenes=True, geneset=None):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    results = []
    for gene_set_library, label in gene_set_libraries.items():
        response = requests.get(
                    ENRICHR_URL +
                       query_string % (user_list_id, gene_set_library)
                )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)
        resultDataframe = pd.DataFrame(data[gene_set_library], columns=[
                                       'rank', 'term_name', 'pvalue', 'zscore', 'combined_score', 'overlapping_genes', 'FDR', 'old_pvalue', 'old_FDR'])
        selectedColumns = ['term_name', 'zscore', 'combined_score', 'pvalue', 'FDR'] if not overlappingGenes else [
            'term_name', 'zscore', 'combined_score', 'FDR', 'pvalue', 'overlapping_genes']
        resultDataframe = resultDataframe.loc[:, selectedColumns]
        resultDataframe['gene_set_library'] = label
        resultDataframe['log10P'] = -np.log10(resultDataframe['pvalue'])
        results.append(resultDataframe)
    concatenatedDataframe = pd.concat(results)
    if geneset:
        concatenatedDataframe['geneset'] = geneset
    return concatenatedDataframe



def get_enrichr_results_by_library(enrichr_results, signature_label, plot_type='interactive', library_type='go', version='2018', sort_results_by='pvalue'):

    # Libraries
    if library_type == 'go':
        go_version = str(version)
        libraries = {
            'GO_Biological_Process_'+go_version: 'Gene Ontology Biological Process ('+go_version+' version)',
            'MGI_Mammalian_Phenotype_Level_4_2019': 'MGI Mammalian Phenotype Level 4 2019'
        }
    elif library_type == "pathway":
        # Libraries
        libraries = {
            'KEGG_2019_Human': 'KEGG Pathways',
        }
    elif library_type == "celltype":
        # Libraries
        libraries = {
            'ARCHS4_Tissues': 'ARCHS4 Tissues',
            'Human_Gene_Atlas': 'Human Gene Atlas',
            'Descartes_Cell_Types_and_Tissue_2021': 'Descartes Cell Types'
        }
    elif library_type=="disease":
        libraries = {
            'GWAS_Catalog_2019': 'GWAS Catalog',
        }
    # Get Enrichment Results
    enrichment_results = {geneset: get_enrichr_results(enrichr_results[geneset]['userListId'], gene_set_libraries=libraries, geneset=geneset) for geneset in ['upregulated', 'downregulated']}
    enrichment_results['signature_label'] = signature_label
    enrichment_results['plot_type'] = plot_type
    enrichment_results['sort_results_by'] = sort_results_by

    # Return
    return enrichment_results
        
    # Get Enrichment Results
    enrichment_results = {geneset: get_enrichr_results(enrichr_results[geneset]['userListId'], gene_set_libraries=libraries, geneset=geneset) for geneset in ['upregulated', 'downregulated']}
    enrichment_results['signature_label'] = signature_label
    enrichment_results['plot_type'] = plot_type
    enrichment_results['sort_results_by'] = sort_results_by

    # Return
    return enrichment_results

# enrichment analysis for uploaded gmt
def get_library(filename):
    # processes library data
    raw_library_data = []
    library_data = []

    
    with open(filename, "r") as f:
        for line in f.readlines():
            raw_library_data.append(line.split("\t\t"))
    name = []
    gene_list = []

    for i in range(len(raw_library_data)):
        name += [raw_library_data[i][0]]
        raw_genes = raw_library_data[i][1].replace('\t', ' ')
        gene_list += [raw_genes[:-1]]

    library_data = [list(a) for a in zip(name, gene_list)]
    
    return library_data

def library_to_dict(library_data):
    dictionary = dict()
    for i in range(len(library_data)):
        row = library_data[i]
        dictionary[row[0]] = row[1].split(" ")
    return dictionary

def get_library_iter(library_data):
    for term in library_data.keys():
        single_set = library_data[term]
        yield term, single_set

def get_enrichment_results(items, library_data):
    return sorted(enrich_crisp(items, get_library_iter(library_data), n_background_entities=20000, preserve_overlap=True), key=lambda r: r[1].pvalue)

# Call enrichment results and return a plot and dataframe for Scatter Plot
def get_values(obj_list):
    pvals = []
    odds_ratio = []
    n_overlap = []
    overlap = []
    for i in obj_list:
        pvals.append(i.pvalue)
        odds_ratio.append(i.odds_ratio)
        n_overlap.append(i.n_overlap)
        overlap.append(i.overlap)
    return pvals, odds_ratio, n_overlap, overlap

def get_qvalue(p_vals):
    r = multipletests(p_vals, method="fdr_bh")
    return r[1]


def enrichment_analysis(items, library_data):    
    all_results = get_enrichment_results(items, library_data)
    unzipped_results = list(zip(*all_results))
    pvals, odds_ratio, n_overlap, overlap = get_values(unzipped_results[1])
    df = pd.DataFrame({"term_name":unzipped_results[0], "p value": pvals, \
                       "odds_ratio": odds_ratio, "n_overlap": n_overlap, "overlap": overlap})
    df["-log(p value)"] = -np.log10(df["p value"])
    df["q value"] = get_qvalue(df["p value"].tolist())
    return [list(unzipped_results[0])], [pvals], df

def get_enrichr_result_tables_by_library(enrichr_results, signature_label, library_type='tf'):

    # Libraries
    if library_type == 'tf':
        # Libraries
        libraries = {
            'ChEA_2016': 'ChEA 2016',
            'ENCODE_TF_ChIP-seq_2015': 'ENCODE TF',
            'ARCHS4_TFs_Coexp': 'ARCHS4 coexpressed TF'
        }
    elif library_type == "ke":
        # Libraries
        libraries = {
            'KEA_2015': 'KEA 2015',
            'ARCHS4_Kinases_Coexp': 'ARCHS4 coexpressed kinases'
        }
    elif library_type == "mirna":
        libraries = {
        'TargetScan_microRNA_2017': 'TargetScan',
        'miRTarBase_2017': 'miRTarBase'
        }


    # Initialize results
    results = []

    # Loop through genesets
    for geneset in ['upregulated', 'downregulated']:

        # Append ChEA results
        enrichment_dataframe = get_enrichr_results(enrichr_results[geneset]['userListId'], gene_set_libraries=libraries, geneset=geneset)
        results.append(enrichment_dataframe)

    # Concatenate results
    enrichment_dataframe = pd.concat(results)

    return {'enrichment_dataframe': enrichment_dataframe, 'signature_label': signature_label}
    

def plot_library_barchart(enrichr_results, gene_set_library, signature_label, case_name, sort_results_by='pvalue', nr_genesets=15, height=400, plot_type='interactive'):
    sort_results_by = 'log10P' if sort_results_by == 'pvalue' else 'combined_score'
    fig = tools.make_subplots(rows=1, cols=2, print_grid=False)
    for i, geneset in enumerate(['upregulated', 'downregulated']):
        # Get dataframe
        enrichment_dataframe = enrichr_results[geneset]
        plot_dataframe = enrichment_dataframe[enrichment_dataframe['gene_set_library'] == gene_set_library].sort_values(sort_results_by, ascending=False).iloc[:nr_genesets].iloc[::-1]

        # Format
        n = 7
        plot_dataframe['nr_genes'] = [len(genes) for genes in plot_dataframe['overlapping_genes']]
        plot_dataframe['overlapping_genes'] = ['<br>'.join([', '.join(genes[i:i+n]) for i in range(0, len(genes), n)]) for genes in plot_dataframe['overlapping_genes']]

        # Get Bar
        bar = go.Bar(
            x=plot_dataframe[sort_results_by],
            y=plot_dataframe['term_name'],
            orientation='h',
            name=geneset.title(),
            showlegend=False,
            hovertext=['<b>{term_name}</b><br><b>P-value</b>: <i>{pvalue:.2}</i><br><b>FDR</b>: <i>{FDR:.2}</i><br><b>Z-score</b>: <i>{zscore:.3}</i><br><b>Combined score</b>: <i>{combined_score:.3}</i><br><b>{nr_genes} Genes</b>: <i>{overlapping_genes}</i><br>'.format(**rowData) for index, rowData in plot_dataframe.iterrows()],
            hoverinfo='text',
            marker={'color': '#FA8072' if geneset == 'upregulated' else '#87CEFA'}
        )
        fig.append_trace(bar, 1, i+1)

        # Get text
        text = go.Scatter(
            x=[max(bar['x'])/50 for x in range(len(bar['y']))],
            y=bar['y'],
            mode='text',
            hoverinfo='none',
            showlegend=False,
            text=['*<b>{}</b>'.format(rowData['term_name']) if rowData['FDR'] < 0.1 else '{}'.format(
                rowData['term_name']) for index, rowData in plot_dataframe.iterrows()],
            textposition="middle right",
            textfont={'color': 'black'}
        )
        fig.append_trace(text, 1, i+1)

    # Get annotations
    annotations = [
        {'x': 0.25, 'y': 1.06, 'text': '<span style="color: #FA8072; font-size: 10pt; font-weight: 600;">Up-regulated in ' +
            case_name+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},
        {'x': 0.75, 'y': 1.06, 'text': '<span style="color: #87CEFA; font-size: 10pt; font-weight: 600;">Down-regulated in ' +
            case_name+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}
    ] if signature_label else []

    # Get title
    title = signature_label + ' | ' + gene_set_library

    fig['layout'].update(height=height, title='<b>{}</b>'.format(title),
                         hovermode='closest', annotations=annotations)
    fig['layout']['xaxis1'].update(domain=[0, 0.49], title='-log10P' if sort_results_by == 'log10P' else 'Enrichment score')
    fig['layout']['xaxis2'].update(domain=[0.51, 1], title='-log10P' if sort_results_by == 'log10P' else 'Enrichment score')
    fig['layout']['yaxis1'].update(showticklabels=False)
    fig['layout']['yaxis2'].update(showticklabels=False)
    fig['layout']['margin'].update(l=0, t=65, r=0, b=35)
    fig.show()



def results_table(enrichment_dataframe, source_label, target_label, table_counter):

    # Get libraries
    for gene_set_library in enrichment_dataframe['gene_set_library'].unique():

        # Get subset
        enrichment_dataframe_subset = enrichment_dataframe[enrichment_dataframe['gene_set_library'] == gene_set_library].copy()

        # Get unique values from source column
        enrichment_dataframe_subset[source_label] = [x.split('_')[0] for x in enrichment_dataframe_subset['term_name']]
        enrichment_dataframe_subset = enrichment_dataframe_subset.sort_values(['FDR', 'pvalue']).rename(columns={'pvalue': 'P-value'}).drop_duplicates(source_label)

        # Add links and bold for significant results
        # if " " in enrichment_dataframe_subset[source_label][0]:
        enrichment_dataframe_subset[source_label] = ['<a href="http://www.mirbase.org/cgi-bin/query.pl?terms={}" target="_blank">{}</a>'.format(x.split(" ")[0], x) if '-miR-' in x else '<a href="http://maayanlab.cloud/Harmonizome/gene/{}" target="_blank">{}</a>'.format(x.split(" ")[0], x)for x in enrichment_dataframe_subset[source_label]]
          
        # else:
        #     enrichment_dataframe_subset[source_label] = ['<a href="http://www.mirbase.org/cgi-bin/query.pl?terms={x}" target="_blank">{x}</a>'.format(**locals()) if '-miR-' in x else '<a href="http://maayanlab.cloud/Harmonizome/gene/{x}" target="_blank">{x}</a>'.format(**locals())for x in enrichment_dataframe_subset[source_label]]
        enrichment_dataframe_subset[source_label] = [rowData[source_label].replace('target="_blank">', 'target="_blank"><b>').replace('</a>', '*</b></a>') if rowData['FDR'] < 0.05 else rowData[source_label] for index, rowData in enrichment_dataframe_subset.iterrows()]

        # Add rank
        enrichment_dataframe_subset['Rank'] = ['<b>'+str(x+1)+'</b>' for x in range(len(enrichment_dataframe_subset.index))]

        # Add overlapping genes with tooltip
        enrichment_dataframe_subset['nr_overlapping_genes'] = [len(x) for x in enrichment_dataframe_subset['overlapping_genes']]
        enrichment_dataframe_subset['overlapping_genes'] = [', '.join(x) for x in enrichment_dataframe_subset['overlapping_genes']]
        enrichment_dataframe_subset[target_label.title()] = ['{nr_overlapping_genes} {geneset} '.format(**rowData)+target_label+'s' for index, rowData in enrichment_dataframe_subset.iterrows()]
        # enrichment_dataframe[target_label.title()] = ['<span class="gene-tooltip">{nr_overlapping_genes} {geneset} '.format(**rowData)+target_label+'s<div class="gene-tooltip-text">{overlapping_genes}</div></span>'.format(**rowData) for index, rowData in enrichment_dataframe.iterrows()]

        # Convert to HTML
        pd.set_option('max.colwidth', -1)
        html_table = enrichment_dataframe_subset.head(50)[['Rank', source_label, 'P-value', 'FDR', target_label.title()]].to_html(escape=False, index=False, classes='w-100')
        html_results = '<div style="max-height: 200px; overflow-y: scroll;">{}</div>'.format(html_table)

        # Add CSS
        display(HTML('<style>.w-100{width: 100%;} .text-left th{text-align: left !important;}</style>'))
        display(HTML('<style>.slick-cell{overflow: visible;}.gene-tooltip{text-decoration: underline; text-decoration-style: dotted;}.gene-tooltip .gene-tooltip-text{visibility: hidden; position: absolute; left: 60%; width: 250px; z-index: 1000; text-align: center; background-color: black; color: white; padding: 5px 10px; border-radius: 5px;} .gene-tooltip:hover .gene-tooltip-text{visibility: visible;} .gene-tooltip .gene-tooltip-text::after {content: " ";position: absolute;bottom: 100%;left: 50%;margin-left: -5px;border-width: 5px;border-style: solid;border-color: transparent transparent black transparent;}</style>'))

        # Display table
        display(HTML(html_results))
        # Display gene set
        # display_object(table_counter, gene_set_library, istable=True)
        if source_label == "Transcription Factor":
            additional_description = ". The table contains scrollable tables displaying the results of the Transcription Factor (TF) enrichment analysis generated using Enrichr. Every row represents a TF; significant TFs are highlighted in bold."
        elif source_label == "Kinase":
            additional_description = ". The table contains browsable tables displaying the results of the Protein Kinase (PK) enrichment analysis generated using Enrichr. Every row represents a PK; significant PKs are highlighted in bold."    
        elif source_label == "miRNA":
            additional_description = ". The figure contains browsable tables displaying the results of the miRNA enrichment analysis generated using Enrichr. Every row represents a miRNA; significant miRNAs are highlighted in bold."
        display_object(table_counter, gene_set_library+additional_description, istable=True)
        display(create_download_link(enrichment_dataframe_subset, filename="Enrichment_analysis_{}_{}.csv".format(source_label, gene_set_library)))
        table_counter += 1
        
    return table_counter

def display_table(analysis_results, source_label, table_counter):
    
    # Plot Table
    return results_table(analysis_results['enrichment_dataframe'].copy(), source_label=source_label, target_label='target', table_counter=table_counter)
