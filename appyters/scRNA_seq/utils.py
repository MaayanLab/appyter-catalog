# Basic libraries
import pandas as pd
import requests, json
import time
import numpy as np
import warnings
import re
import random
from collections import defaultdict

# Visualization
import scipy.stats as ss
import plotly
from plotly import tools
import plotly.express as px
import plotly.graph_objs as go
import matplotlib
import matplotlib.pyplot as plt; plt.rcdefaults()
from matplotlib import rcParams
from matplotlib.lines import Line2D
from matplotlib_venn import venn2, venn3
import IPython
from IPython.display import HTML, display, Markdown, IFrame, FileLink
from itertools import combinations, permutations
from scipy import stats
import seaborn as sns
# Data analysis
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform
from sklearn import cluster
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
import umap
from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
from magic import MAGIC as MG
import scanpy as sc
import anndata
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
from maayanlab_bioinformatics.dge.limma_voom import limma_voom_differential_expression
import pandas as pd
import sys, h5py, time
import scanpy as sc
import anndata
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
import umap.umap_ as umap
from sklearn.decomposition import NMF
from statsmodels.stats.multitest import multipletests

from IPython.display import display, HTML
from maayanlab_bioinformatics.enrichment.crisp import enrich_crisp, fisher_overlap

# Bokeh
from bokeh.io import output_notebook
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, CustomJS, ColumnDataSource, Span, Select, Legend, PreText, Paragraph, LinearColorMapper, ColorBar, CategoricalColorMapper
from bokeh.layouts import layout, row, column, gridplot
from bokeh.palettes import all_palettes
import colorcet as cc
from bokeh.palettes import Category20

from plotly.offline import init_notebook_mode
init_notebook_mode(connected = False)
output_notebook()


pd.set_option('display.max_columns', 1000)  
pd.set_option('display.max_rows', 1000)


def check_files(fname):
    if fname == "":
        raise IOError
    if fname.endswith(".txt") == False and fname.endswith(".csv") ==False and fname.endswith(".tsv")==False:
        raise IOError
def check_df(df, col):
    if col not in df.columns:
        raise IOError

def load_seurat_files(mtx_filename, gene_filename, barcodes_filename):
    
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
    adata.var['gene_symbols'] = genes.iloc[:, 1].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")
    
    
    adata.obs['barcode'] = cells
    adata.obs_names = cells
    adata.obs_names_make_unique(join="-")
    return adata

def load_metadata(adata, meta_data_filename, meta_class_column_name):
    if meta_data_filename != "":
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

    else:
        meta_class_column_name = "Class"
        adata.obs[meta_class_column_name] = ["Class0"]*adata.n_obs
        adata.var_names_make_unique()
    
    return adata, meta_class_column_name

def load_data(dataset_name, rnaseq_data_filename, mtx_data_filename, gene_data_filename, barcode_data_filename, meta_data_filename=None, meta_class_column_name=None, table_counter=1):
    adata = None
    if rnaseq_data_filename != "":
        check_files(rnaseq_data_filename)
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
    elif mtx_data_filename != "":
        adata = load_seurat_files(mtx_data_filename, gene_data_filename, barcode_data_filename)
    if adata is not None:
        # load meta data
        adata, meta_class_column_name = load_metadata(adata, meta_data_filename, meta_class_column_name)    
    
        # add batch info
        if meta_class_column_name == "":
            meta_class_column_name = "batch"
        else:
            adata.obs = adata.obs.rename(columns={meta_class_column_name: 'class'})
            meta_class_column_name = "class"
            
        adata.obs["batch"] = dataset_name
        adata.obs.index = adata.obs.index + "-" + dataset_name
        
        display_statistics(adata, "### Statistics of data ###") 
    return adata, table_counter, meta_class_column_name
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

def display_object(counter, caption, df=None, istable=True, subcounter=""):
    if df is not None:
        display(df)
    if istable == True:
        display(Markdown("*Table {}{}. {}*".format(counter, subcounter, caption)))
    else:
        display(Markdown("*Figure {}{}. {}*".format(counter, subcounter, caption)))
    counter += 1
    return counter


def display_statistics(data, description=""):
    print(description)
    print("Sample size:", data.n_obs)
    print("Feature size:", data.n_vars)
    
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





#############################################
########## 2. Plot
#############################################


def normalize(adata, normalization_method, log_normalization):
    tmp_adata = adata.copy()
    if normalization_method == "Seurat":
        sc.pp.filter_cells(tmp_adata, min_genes=200)
        sc.pp.filter_genes(tmp_adata, min_cells=3)
        sc.pp.normalize_total(tmp_adata, target_sum=1e4)        
        if log_normalization:
            sc.pp.log1p(tmp_adata)
        sc.pp.scale(tmp_adata, max_value=10)
    elif normalization_method == "Zheng17":
        sc.pp.recipe_zheng17(adata_merged, log=log_normalization, plot=False)
    elif normalization_method == "Weinreb17":
        sc.pp.recipe_weinreb17(adata_merged, log=log_normalization)
    return tmp_adata
    


def run_magic(dataset, solver='exact'):
    # Run imputation
    dataset.uns['magic'] = normalize_magic(dataset.to_df(), solver=solver).T
    return dataset
def normalize_magic(dataset, solver='exact'):
    
    magic_op = MG(solver=solver)
    data_magic = magic_op.fit_transform(dataset)
    return data_magic.transpose()
    

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
    clustergrammer_url = 'https://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
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

    
    if cluster == True:
        # cluster 0 vs rest
        sc.tl.rank_genes_groups(dataset, meta_class_column_name, method='t-test', use_raw=True)
            
        for cls1 in classes:
            signature_label = f"{cls1} vs. rest"
            print("Analyzing.. {} using {}".format(signature_label, method))
            cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls1, :].index.tolist() #case
            non_cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]!=cls1, :].index.tolist() #control
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
                signature = signature.sort_values("logFC", ascending=False)
            elif method == "DESeq2":
                DESeq2 = robjects.r['deseq2']
                DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(non_cls1_sample_ids), pandas2ri.conversion.py2rpy(cls1_sample_ids)))

                signature = pd.DataFrame(DESeq2_results[0])
                signature.index = DESeq2_results[1]
                signature = signature.sort_values("log2FoldChange", ascending=False)
            elif method == "wilcoxon":   
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls1).set_index('names').sort_values('pvals', ascending=True)
                dedf = dedf.replace([np.inf, -np.inf], np.nan).dropna()              
                dedf = dedf.sort_values("logfoldchanges", ascending=False)
                signature = dedf
                
            signatures[signature_label] = signature
    else:
        sc.tl.rank_genes_groups(dataset, meta_class_column_name, method='wilcoxon', use_raw=True)
                
        for cls1, cls2 in permutations(classes, 2):
            signature_label = " vs. ".join([cls1, cls2])
            print("Analyzing.. {} using {}".format(signature_label, method))
            cls1_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls1, :].index.tolist() #control
            cls2_sample_ids = meta_df.loc[meta_df[meta_class_column_name]==cls2, :].index.tolist() #case
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
                signature = signature.sort_values("logFC", ascending=False)
            elif method == "DESeq2":
                DESeq2 = robjects.r['deseq2']
                DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))

                signature = pd.DataFrame(DESeq2_results[0])
                signature.index = DESeq2_results[1]
                signature = signature.sort_values("log2FoldChange", ascending=False)
            elif method == "wilcoxon":   
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls2).set_index('names').sort_values('pvals', ascending=True)
                dedf = dedf.replace([np.inf, -np.inf], np.nan).dropna()
                dedf = dedf.sort_values("logfoldchanges", ascending=False)
                signature = dedf
            signatures[signature_label] = signature
    return signatures




def submit_enrichr_geneset(geneset, label):
    ENRICHR_URL = 'https://amp.pharm.mssm.edu/Enrichr/addList'
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


def run_enrichr(signature, signature_label, geneset_size=500, fc_colname = 'logFC'):

    # Sort signature
    up_signature = signature[signature[fc_colname] > 0]
    down_signature = signature[signature[fc_colname] < 0]
    
    # Get genesets
    genesets = {
        'upregulated': up_signature.index[:geneset_size],
        'downregulated': down_signature.index[:geneset_size:]
    }

    # Submit to Enrichr
    enrichr_ids = {geneset_label: submit_enrichr_geneset(geneset=geneset, label=signature_label+', '+geneset_label+', from scRNA-seq Appyter') for geneset_label, geneset in genesets.items()}
    enrichr_ids['signature_label'] = signature_label
    return enrichr_ids

def get_enrichr_results(user_list_id, gene_set_libraries, overlappingGenes=True, geneset=None):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
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
            'HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression': 'HuBMAP ASCT+B Cell Type'
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


def get_enrichr_result_tables_by_library(enrichr_results, signature_label, library_type='tf'):

    # Libraries
    if library_type == 'tf':
        # Libraries
        libraries = {
            'ChEA_2016': 'ChEA (experimentally validated targets)',
        }
    elif library_type == "ke":
        # Libraries
        libraries = {
            'KEA_2015': 'KEA (experimentally validated targets)',
            'ARCHS4_Kinases_Coexp': 'ARCHS4 (coexpressed genes)'
        }
    elif library_type == "mirna":
        libraries = {
        'TargetScan_microRNA_2017': 'TargetScan (experimentally validated targets)',
        'miRTarBase_2017': 'miRTarBase (experimentally validated targets)'
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
        dictionary[row[0]] = [x.upper() for x in row[1].split(" ")]
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
    items = [x.upper() for x in items]
    all_results = get_enrichment_results(items, library_data)
    unzipped_results = list(zip(*all_results))
    if len(unzipped_results)>0:
        pvals, odds_ratio, n_overlap, overlap = get_values(unzipped_results[1])
        df = pd.DataFrame({"term_name":unzipped_results[0], "pvalue": pvals, \
                           "odds_ratio": odds_ratio, "n_overlap": n_overlap, "overlap": overlap})
        df["-log(p value)"] = -np.log10(df["pvalue"])
        df["q value"] = get_qvalue(df["pvalue"].tolist())
        return [list(unzipped_results[0])], [pvals], df
    else:
        raise IOError()


import hashlib
def str_to_int(string, mod):
    string = re.sub(r"\([^()]*\)", "", string).strip()
    byte_string = bytearray(string, "utf8")
    return int(hashlib.sha256(byte_string).hexdigest(), base=16)%mod

def plot_scatter(umap_df, values_dict, option_list, sample_names, caption_text, category_list_dict=None, location='right', category=True, dropdown=False, figure_counter=0, additional_info=None):
    
    # init plot 
    if additional_info is not None:
        source = ColumnDataSource(data=dict(x=umap_df["x"], y=umap_df["y"], values=values_dict[option_list[0]], 
                                        names=sample_names, info=additional_info[option_list[0]]))
    else:
        source = ColumnDataSource(data=dict(x=umap_df["x"], y=umap_df["y"], values=values_dict[option_list[0]], 
                                        names=sample_names))
    # node size
    if umap_df.shape[0] > 1000:
        node_size = 2
    else:
        node_size = 6
        
    if location == 'right':
        plot = figure(plot_width=1000, plot_height=800)   
    else:
        plot = figure(plot_width=1000, plot_height=1000+20*len(category_list_dict[option_list[0]]))   
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
        scatter = plot.scatter('x', 'y', size=node_size, source=source, color={'field': 'values', 'transform': color_mapper}, legend_field="values")
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
        plot.scatter('x', 'y', size=node_size,  source=source, color={'field': 'values', 'transform': color_mapper})
    
    if additional_info is not None:
            tooltips = [
            ("Sample", "@names"),
            ("Value", "@values"),
            ("p-value", "@info")
        ]
    else:
        tooltips = [
            ("Sample", "@names"),
            ("Value", "@values"),
        ]
    plot.add_tools(HoverTool(tooltips=tooltips))
    plot.output_backend = "webgl"
    
    plot.xaxis.axis_label = "UMAP_1"
    plot.xaxis.axis_label_text_font_size = "12pt"
    plot.yaxis.axis_label = "UMAP_2"
    plot.yaxis.axis_label_text_font_size = "12pt"
    
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
                                              additional_info=additional_info,\
                                              factors_dict=factors_dict,\
                                              colors_dict=colors_dict,\
                                              plot=plot,\
                                              scatter=scatter,
                                              caption_text=caption_text
                                             ), code="""        
                const val = cb_obj.value;                    
                source.data.values = values_dict[val]  
                if (additional_info != null) {
                    source.data.info = additional_info[val]
                }
                color_mapper.factors = category_list_dict[val]
                color_mapper.palette = colors_dict[val]
                plot.legend = unique_category_dict[val]
                pre.text = "Figure "+figure_counter+". "+caption_text+val+"."; 
                plot.height = 1000+20*(category_list_dict[val].length)
                source.change.emit();
            """)
        else:
            callback_adt = CustomJS(args=dict(source=source, \
                                              pre=pre, \
                                              values_dict=values_dict, \
                                              additional_info=additional_info,\
                                              figure_counter=figure_counter,
                                              caption_text=caption_text), code="""        
                const val = cb_obj.value;    
                source.data.values = values_dict[val]
                if (additional_info != null) {
                    source.data.info = additional_info[val]
                }
                
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




def clustering(adata, dataset, umap_n_neighbors, umap_min_dist, bool_plot, figure_counter, batch_correction=True):
    # clustering
    
    if batch_correction== True:
        sc.external.pp.bbknn(adata, batch_key="batch")
    else:
        sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors)
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.umap(adata, min_dist=umap_min_dist)
    
    
    # sort by clusters
    new_order = adata.obs.sort_values(by='leiden').index.tolist()
    adata = adata[new_order, :]
    adata.obs['leiden'] = 'Cluster '+adata.obs['leiden'].astype('object')
    
    if bool_plot == True:
        umap_df = pd.DataFrame(adata.obsm['X_umap'])
        umap_df.columns = ['x', 'y']

        values_dict = dict()
        values_dict["Cluster"] = adata.obs["leiden"].values
        category_list_dict = dict()
        category_list_dict["Cluster"] = list(sorted(adata.obs["leiden"].unique()))
        figure_counter = plot_scatter(umap_df, values_dict, ["Cluster"], adata.obs.index.tolist(), "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False, figure_counter=figure_counter)

        display(create_download_link(adata.obs["leiden"], filename=f"clustering_{dataset}.csv"))
    return adata, figure_counter


def differential_gene_expression_analysis(adata, diff_gex_method, enrichment_groupby, meta_class_column_name, table_counter):
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
        sort_genes_by = "scores"
        ascending = False
        
        
    if enrichment_groupby == "user_defined_class":
        classes = adata.obs[meta_class_column_name].unique().tolist()
        bool_cluster=False
        if len(classes) < 2:
            print("Warning: Please provide at least 2 classes in the metadata")
        
    elif enrichment_groupby == "Cluster":
        meta_class_column_name = "leiden"
        classes = sorted(adata.obs["leiden"].unique().tolist())
        classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
        bool_cluster=True
    else:
        meta_class_column_name = enrichment_groupby
        classes = sorted(adata.obs[meta_class_column_name].unique().tolist())
        classes.sort()
        bool_cluster=True
        
    if len(classes) > 5 and adata.n_obs > 5000:
        if diff_gex_method == "wilcoxon":
            print('Warning: There are too many cells/clusters. It cannot execute the analysis code for the data. The appyter randomly select 5000 samples. If you want to execute it with the whole data, please run it locally.')
        else:
            print('Warning: There are too many cells/clusters. It cannot execute the analysis code for the data. The appyter switched to Wilcoxon rank-sum method and randomly select 5000 samples. If you want to execute it with the whole data, please run it locally.')
            diff_gex_method = "wilcoxon"
        # randomly select 5K samples
        random_selected_samples = random.sample(adata.obs.index.tolist(), 5000)
        adata_random_sampled = adata[random_selected_samples, :]
        signatures = get_signatures(classes, adata_random_sampled, method=diff_gex_method, meta_class_column_name=meta_class_column_name, cluster=bool_cluster)
    else:
        signatures = get_signatures(classes, adata, method=diff_gex_method, meta_class_column_name=meta_class_column_name, cluster=bool_cluster)
    
    if len(classes) > 1:
        for label, signature in signatures.items():
            if bool_cluster == True:
                case_label = label.split(" vs. ")[0]
            else:
                case_label = label.split(" vs. ")[1]
            signature = signature.sort_values(by=sort_genes_by, ascending=ascending)
            table_counter = display_object(table_counter, f"Top 5 Differentially Expressed Genes in {label} (up-regulated in {case_label})", signature.head(5), istable=True)
            display(create_download_link(signature, filename="DEG_{}.csv".format(label)))
    return signatures, bool_cluster, table_counter



def visualize_enrichment_analysis(adata, signatures, meta_class_column_name, diff_gex_method, enrichr_libraries_filename, enrichr_libraries, enrichment_groupby, libraries_tab, gene_topk, bool_cluster, bool_plot, figure_counter, table_counter):

    topk_enrichment_terms = 3
    if diff_gex_method == "characteristic_direction":
        fc_colname = "CD-coefficient"
    elif diff_gex_method == "limma":
        fc_colname = "logFC"
    elif diff_gex_method == "edgeR":
        fc_colname = "logFC"
    elif diff_gex_method == "DESeq2":
        fc_colname = "log2FoldChange"
    elif diff_gex_method == "wilcoxon":
        fc_colname = "logfoldchanges"
        
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

            results['enrichr'][label] = run_enrichr(signature=signature, signature_label=label, fc_colname=fc_colname, geneset_size=gene_topk)
            display(Markdown("*Enrichment Analysis Result: {} (Up-regulated in {})*".format(label, case_name)))
            display_link("https://amp.pharm.mssm.edu/Enrichr/enrich?dataset={}".format(results['enrichr'][label]["upregulated"]["shortId"]))
            display(Markdown("*Enrichment Analysis Result: {} (Down-regulated in {})*".format(label, case_name)))
            display_link("https://amp.pharm.mssm.edu/Enrichr/enrich?dataset={}".format(results['enrichr'][label]["downregulated"]["shortId"]))
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
            up_signature = signature[signature[fc_colname] > 0]
            up_genes = [x.upper() for x in up_signature.index[:gene_topk].tolist()]

            results['user_defined_enrichment'][label] = dict()
            try:
                _, _, results['user_defined_enrichment'][label]['enrichment_dataframe'] = enrichment_analysis(up_genes, user_library)
                results['user_defined_enrichment'][label]['enrichment_dataframe']['gene_set_library'] = enrichr_libraries_filename
            except:
                pass
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
    topk_enriched_terms_dict = defaultdict(list)
    for label, signature in signatures.items():
        if bool_cluster == True:
            cluster_names = [label.split(" vs. ")[0]]
        else:
            cluster_names = [label.split(" vs. ")[1]]
        
        for key in results.keys():
            if key.endswith("enrichment") == False:
                continue
            enrichment_results = results[key][label]
            meta_df = adata.obs
            for cluster_name in cluster_names:
                for direction in ['upregulated']:
                    if direction in enrichment_results:
                        enrichment_dataframe = enrichment_results[direction]
                    elif "enrichment_dataframe" in enrichment_results:
                        enrichment_dataframe = enrichment_results["enrichment_dataframe"]
                    else:
                        continue
                    if enrichment_dataframe.empty == True:
                        raise Exception("Enrichment analysis returns empty results. Please check if your data contains proper gene names.")
                    libraries = enrichment_dataframe['gene_set_library'].unique() 
                    
                    for library in libraries:
                        enrichment_dataframe_library = enrichment_dataframe[enrichment_dataframe['gene_set_library']==library]
                        
                        topk_enriched_terms = enrichment_dataframe_library.iloc[:topk_enrichment_terms]
                        selected_columns = ['term_name', 'pvalue']
                        topk_enriched_terms = topk_enriched_terms[selected_columns]
                        topk_enriched_terms.columns = ["Enriched Term", "pvalue"]
                        topk_enriched_terms.insert(0, 'Rank', [*range(1, len(topk_enriched_terms)+1)])
                        topk_enriched_terms.insert(0, 'Class', cluster_name)
                        
                        topk_enriched_terms_dict[library].append(topk_enriched_terms)
                        
                        top_term = enrichment_dataframe_library.iloc[0]['term_name']
                        if library not in meta_df.columns:
                            meta_df.insert(0, library, np.nan)
                        meta_df[library] = meta_df[library].astype("object")
                        if bool_cluster == True:
                            meta_df.loc[meta_df[col_name]==cluster_name, library] = top_term
                        else:
                            meta_df.loc[meta_df[meta_class_column_name]==cluster_name, library] = top_term
                        library_option_list.add(library)
                        
                        
                        for topk_i in range(topk_enrichment_terms):
                            if topk_i < topk_enriched_terms.shape[0]: # in case # of enrichment terms is less than topk_enrichment_terms
                                new_col_name = f"{library}_{topk_i}"
                                top_term = topk_enriched_terms.iloc[topk_i]['Enriched Term']

                                if new_col_name not in meta_df.columns:
                                    meta_df.insert(0, new_col_name, np.nan)
                                meta_df[new_col_name] = meta_df[new_col_name].astype("object")
                                meta_df.loc[meta_df[col_name]==cluster_name, new_col_name] = top_term

    
    library_option_list = list(library_option_list)

    for library in library_option_list: 
        top_enriched_term_df = pd.concat(topk_enriched_terms_dict[library])
        top_enriched_term_df = top_enriched_term_df.set_index(["Class", "Rank"])
        
        if libraries_tab == 'Yes':
            table_counter = display_object(table_counter, f"Top 3 Enriched Terms for each cluster/class in library {library}. To see more results, please use Enrichr links above.", df=top_enriched_term_df, istable=True)
        else:
            table_counter = display_object(table_counter, f"Top 3 Enriched Terms for each cluster/class in library {library}.", df=top_enriched_term_df, istable=True)
        display(create_download_link(top_enriched_term_df, filename=f"Top{topk_enrichment_terms}_Enriched_Terms_{library}.csv"))
    
    # umap info into dataframe 
    umap_df = pd.DataFrame(adata.obsm['X_umap'])
    umap_df.columns = ['x', 'y']

    option_list = library_option_list  
    if enrichment_groupby != "batch":
        option_list.append("leiden")
    else:
        option_list.append("batch")
    adata_norm_selected = adata.obs[option_list].astype("object").fillna("NaN")
    
    values_dict = dict(zip(adata_norm_selected.T.index.tolist(), adata_norm_selected.T.values))
    category_list_dict = dict()
    for option in option_list:
        category_list_dict[option] = list(sorted(adata_norm_selected[option].unique()))
    if bool_plot == True:
        figure_counter = plot_scatter(umap_df, values_dict, option_list, adata.obs.index.tolist(), "Scatter plot of the samples. Each dot represents a sample and it is colored by enriched terms in library ", location='below', category_list_dict=category_list_dict, category=True, dropdown=True, figure_counter=figure_counter)
    return adata, option_list, figure_counter, table_counter



from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform

def extract_library_name(path):
    return path.replace("libraries/","").replace(".txt","").replace(".gmt", "").replace("_0","")

def plot_protrack(libraries, tmp_adata, figure_counter=0, fig_size=(1800, 800)):
    
    frequent_genes = get_frequent_genes(tmp_adata.uns['signature'])
    frequent_genes = frequent_genes.intersection(tmp_adata.var.index)
    
    # init figure
    trace_idx = 0
    row_length = len(libraries)+2
    
    row_heights = [0.1]
    row_heights.extend([(1-row_heights[0])/(row_length-1)]*(row_length-1))
    fig = make_subplots(rows=row_length, cols=1, shared_xaxes=True, vertical_spacing=0.00, row_heights=row_heights,
                       )

    # plot dendrogram
    fig, dendrogram_fig, trace_idx, dendrogram_sample_index_list = plot_dendrogram(fig, tmp_adata, frequent_genes, trace_idx)
    
    # get color mapping dictionary
    color_map_dict = get_color_map(libraries, tmp_adata)
        
    i = 2

    for library_filename in libraries:
        # init rug plots
        adata_heatmap = tmp_adata.copy()
        new_library_filename = extract_library_name(library_filename)
        library_name = library_filename+"_0"
        clean_library_filename = {library_name: new_library_filename}
        new_cols = list()
        for col in adata_heatmap.obs.columns:
            if col in clean_library_filename:
                new_cols.append(clean_library_filename[col])
            else:
                new_cols.append(col)
        adata_heatmap.obs.columns = new_cols

        ####### rug plots to show predicted cell types/enriched terms    ########
        fig, trace_idx, i = plot_rug(fig, dendrogram_fig, dendrogram_sample_index_list, adata_heatmap, new_library_filename, color_map_dict, trace_idx, i)
        
        
    # Edit Layout
    fig.update_yaxes(showgrid=False, showline=False, zeroline=False, showticklabels=False, ticks="", row=1, col=1)
    fig.update_layout({'width':fig_size[0], 'height':fig_size[1],
                       'legend_traceorder': "grouped",
                             })
    fig.update_layout(paper_bgcolor='#FFFFFF', plot_bgcolor='#FFFFFF') 

    
    # Plot!
    fig.show() 
    figure_counter = display_object(figure_counter, "Protrack Plot. The plot shows the top predicted cell types for each cluster for each library. The colors are based on Cell Ontology IDs.", istable=False)        
    return figure_counter
    
def plot_rug(fig, dendrogram_fig, dendrogram_sample_index_list, adata_heatmap, new_library_filename, color_map_dict, trace_idx, i):
    rug_df = pd.DataFrame(dendrogram_fig['layout']['xaxis']['tickvals'], index=dendrogram_sample_index_list)
    rug_df.columns = ['x']
    rug_df['y'] = new_library_filename
    rug_df['sample'] = rug_df.index
    rug_df['Prediction'] = adata_heatmap.obs.loc[dendrogram_sample_index_list, new_library_filename].map(str)
    cell_types = list(map(str, rug_df['Prediction'].tolist()))

    subfig = px.scatter(rug_df, x='x', y='y', color='Prediction', hover_name='sample',
                        color_discrete_map=color_map_dict, hover_data=['Prediction','y'])
    subfig.update_layout({'width':800, 'height':30, 'showlegend': False})
    subfig.update_traces(marker=dict(size=12, symbol='line-ns-open',
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))
    for sub in subfig['data']:
        fig.add_trace(sub, row=i, col=1)
        fig['data'][trace_idx]['showlegend']=True
        fig['data'][trace_idx]['legendgroup']=new_library_filename
        fig['data'][trace_idx]['legendgrouptitle_text']=new_library_filename

        trace_idx += 1
    # Update yaxis properties
    fig.update_yaxes(showgrid=False, showline=False, zeroline=False, showticklabels=True, ticks="", row=i, col=1)
    fig.update_xaxes(showgrid=False, showline=False, zeroline=False, showticklabels=False, ticks="", row=i, col=1)
    i += 1
    
    return fig, trace_idx, i
def plot_dendrogram(fig, tmp_adata, frequent_genes, trace_idx=0):
    labels = tmp_adata.obs.index # ordered in cluster
    data_array = tmp_adata.to_df().loc[:,frequent_genes]
    genes = data_array.columns

    # Initialize figure by creating upper dendrogram for samples
    dendrogram_fig = ff.create_dendrogram(data_array, orientation='bottom', labels=list(labels), colorscale=["#000000"], color_threshold=0.0) # fig ordered in dendrogram
    for i in range(len(fig['data'])):
        dendrogram_fig['data'][i]['yaxis'] = 'y2'
    dendrogram_fig.update_layout({'width':1000, 'height':800, 'showlegend':False, 'hovermode': 'closest',})
    dendrogram_fig.update_yaxes(mirror=False, showgrid=False, showline=False, zeroline=False, showticklabels=False, ticks="")
    dendrogram_fig.update_xaxes(mirror=False, showgrid=False, showline=False, zeroline=False, showticklabels=False, ticks="")
    dendrogram_sample_index_list = dendrogram_fig['layout']['xaxis']['ticktext'] # sample names
    for subfig in dendrogram_fig['data']:
        fig.add_trace(subfig, row=1, col=1)
        fig['data'][trace_idx]['showlegend']=False
        trace_idx += 1
    return fig, dendrogram_fig, trace_idx, dendrogram_sample_index_list

def get_frequent_genes(signatures, topk_genes_heatmap=1000):
    # frequent genes for heatmap
    frequent_upregulated_genes = dict()
    
    rank_list = list()
    top1000_rank_list = list()
    for label, signature in signatures.items():
        rank_list.append(signature.rank(ascending=False)["scores"])
        top1000_rank_list.append(signature.rank(ascending=False)["scores"]<topk_genes_heatmap)

    rank_df = pd.concat(rank_list, axis=1)
    rank_df["median"] = rank_df.median(axis=1)

    top1000_rank_df = pd.concat(top1000_rank_list, axis=1)
    top1000_rank_df["sum"] = top1000_rank_df.sum(axis=1)\

    top1000_rank_df["median"] = rank_df["median"]
    top1000_rank_df = top1000_rank_df.sort_values(["sum", "median"], ascending=[False, True])

    sorted_dict = set(top1000_rank_df.iloc[:topk_genes_heatmap,:].index)
    
    return sorted_dict
def get_color_map(library_filenames, tmp_adata):
    # legend color map
    color_list = list()
    color_list.extend(px.colors.qualitative.Plotly)
    color_list.extend(px.colors.qualitative.Alphabet)

    color_map_dict = dict() # pred cell ontology name+ID :color
    color_map_id_dict = dict() # pred cell ontology ID : color

    i = 0
    for library_filename in library_filenames:
        library_name = library_filename+"_0" # top1 prediction

        prediction = tmp_adata.obs[library_name].unique()
        for pred in prediction:
            pred = str(pred)
            cl_id_in_pred = pred.split(":")[-1]
            if cl_id_in_pred == "": #if CL ID is None, use cell ontology name
                cl_id_in_pred = pred

            if cl_id_in_pred not in color_map_id_dict:
                color_map_id_dict[cl_id_in_pred] = color_list[i]
                color_map_dict[pred] = color_list[i]
                i += 1
            else:
                color_map_dict[pred] = color_map_id_dict[cl_id_in_pred]


            if i == len(color_list):
                i = i-len(color_list)
    return color_map_dict
def summary(adata, option_list, table_counter):
    for col in option_list:
        counts = adata.obs[[col]].reset_index().groupby(col).count()
        counts.columns = ['# of Samples']
        counts["Percentage (%)"] = counts['# of Samples']/counts['# of Samples'].sum() * 100
        counts = counts.sort_values("Percentage (%)", ascending=False)
        table_counter = display_object(table_counter, "The number of samples for each category in {}".format(col), counts, istable=True)
    return table_counter
def trajectory_inference(adata, trajectory_method, figure_counter=0):
    node_size = min(100, 120000 / len(adata.obs.index))
    if trajectory_method == "palantir":
        # Run analysis
        pr_res, umap = run_palantir(dataset=adata)
        
        # Display results
        sc.pl.umap(adata, color=['leiden'], size=node_size)
        plot_palantir(pr_res, umap)

    elif trajectory_method == "dpt":
        adata.uns['iroot'] = 0
        sc.pl.umap(adata, color=['leiden'], size=node_size)
        sc.tl.dpt(adata)
        sc.pl.umap(adata, color=['dpt_pseudotime'], size=node_size)
        display_link("draw_graph_fa.pdf", "Download figure")
        
    figure_counter = display_object(figure_counter, "Trajectory inference result using {}. Each point represents an RNA-seq sample. Sample colors are based on pseudotime.".format(trajectory_method), istable=False)
    return adata, figure_counter

def run_palantir(dataset):
    sc.external.tl.phenograph(dataset, clustering_algo="leiden", k=30)
    sc.external.tl.palantir(dataset)
    umap = pd.DataFrame(dataset.obsm['X_umap'], index=dataset.obs_names, columns=['x', 'y'])
    
    # early cell is set based on CD34 expression
    df = dataset.to_df()
    df.columns = df.columns.map(str.upper)
    start_cell = df['CD34'].sort_values(ascending=False).index[0]
    pr_res = sc.external.tl.palantir_results(
        dataset,
        early_cell=start_cell,
        ms_data='X_palantir_multiscale',
        num_waypoints=500,
    )
    return pr_res, umap
def plot_palantir(pr_res, umap):
    # Set up figure
    n_branches = pr_res.branch_probs.shape[1]
    n_cols = 6
    n_rows = int(np.ceil(n_branches / n_cols))
    fig = plt.figure(figsize=[2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(
        n_rows + 2, n_cols, height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows))
    )
    cmap = plt.cm.plasma
    # Pseudotime
    ax = plt.subplot(gs[0:2, 1:3])
    c = pr_res.pseudotime[umap.index]
    ax.scatter(umap.loc[:, "x"], umap.loc[:, "y"], s=3, cmap=plt.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Pseudotime")

    # Entropy
    ax = plt.subplot(gs[0:2, 3:5])
    c = pr_res.entropy[umap.index]
    ax.scatter(umap.loc[:, "x"], umap.loc[:, "y"], s=3, cmap=plt.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Differentiation potential")

def time_series_trajectory_inference(adata, timepoint_labels_column_name, timepoint_labels, figure_counter):
    run_tempora(adata, timepoint_labels_column_name, timepoint_labels)
    display(Image.open("Tempora_plot.jpg"))
    figure_counter = display_object(figure_counter, "Time-series trajectory inference result using Tempora. Tempora visualizes the result as a network, with the piechart at each node representing the composition of cells collected at different time points in the experiment and the arrow connecting each pair of nodes representing lineage relationship between them.", istable=False)
    display(FileLink("Tempora_plot.jpg", result_html_prefix="Download figure"))
    return figure_counter

def results_table(enrichment_dataframe, source_label, target_label, label, table_counter):

    # Get libraries
    for gene_set_library in enrichment_dataframe['gene_set_library'].unique():

        # Get subset
        enrichment_dataframe_subset = enrichment_dataframe[enrichment_dataframe['gene_set_library'] == gene_set_library].copy()

        # Get unique values from source column
        enrichment_dataframe_subset[source_label] = [x.split('_')[0] for x in enrichment_dataframe_subset['term_name']]
        enrichment_dataframe_subset = enrichment_dataframe_subset.sort_values(['FDR', 'pvalue']).rename(columns={'pvalue': 'P-value'}).drop_duplicates(source_label)

        # Add links and bold for significant results
        enrichment_dataframe_subset[source_label] = ['<a href="http://www.mirbase.org/cgi-bin/query.pl?terms={}" target="_blank">{}</a>'.format(x.split(" ")[0], x) if '-miR-' in x else '<a href="http://amp.pharm.mssm.edu/Harmonizome/gene/{}" target="_blank">{}</a>'.format(x.split(" ")[0], x)for x in enrichment_dataframe_subset[source_label]]
          
        # else:
        enrichment_dataframe_subset[source_label] = [rowData[source_label].replace('target="_blank">', 'target="_blank"><b>').replace('</a>', '*</b></a>') if rowData['FDR'] < 0.05 else rowData[source_label] for index, rowData in enrichment_dataframe_subset.iterrows()]

        # Add rank
        enrichment_dataframe_subset['Rank'] = ['<b>'+str(x+1)+'</b>' for x in range(len(enrichment_dataframe_subset.index))]

        # Add overlapping genes with tooltip
        enrichment_dataframe_subset['nr_overlapping_genes'] = [len(x) for x in enrichment_dataframe_subset['overlapping_genes']]
        enrichment_dataframe_subset['overlapping_genes'] = [', '.join(x) for x in enrichment_dataframe_subset['overlapping_genes']]
        enrichment_dataframe_subset[target_label.title()] = ['{nr_overlapping_genes} {geneset} '.format(**rowData)+target_label+'s' for index, rowData in enrichment_dataframe_subset.iterrows()]
        
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
        
        if source_label == "Transcription Factor":
            additional_description = "Enrichment Analysis Results for {} in {}. The table contains scrollable tables displaying the results of the Transcription Factor (TF) enrichment analysis generated using Enrichr. Every row represents a TF; significant TFs are highlighted in bold."
        elif source_label == "Kinase":
            additional_description = "Enrichment Analysis Results for {} in {}. The table contains browsable tables displaying the results of the Protein Kinase (PK) enrichment analysis generated using Enrichr. Every row represents a PK; significant PKs are highlighted in bold."    
        elif source_label == "miRNA":
            additional_description = "Enrichment Analysis Results for {} in {}. The figure contains browsable tables displaying the results of the miRNA enrichment analysis generated using Enrichr. Every row represents a miRNA; significant miRNAs are highlighted in bold."
        display_object(table_counter, additional_description.format(label, gene_set_library), istable=True)
        display(create_download_link(enrichment_dataframe_subset, filename="Enrichment_analysis_{}_{}.csv".format(source_label, gene_set_library)))
        table_counter += 1
        
    return table_counter

def display_table(analysis_results, source_label, label, table_counter):
    
    # Plot Table
    return results_table(analysis_results['enrichment_dataframe'].copy(), source_label=source_label, target_label='target', label=label, table_counter=table_counter)
def CPM(data):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)

    return data


def run_monocle(dataset, color_by='Pseudotime', ordering='de', plot_type='interactive'):
    robjects.r('''
    suppressMessages(library(dplyr))
    suppressMessages(library(monocle))
    suppressMessages(library(tibble))
    suppressMessages(require(Matrix))
    suppressMessages(require(VGAM))
    suppressMessages(require(igraph))

    # Make a CellDataSet object
    # @expr_df: CPM expression data.frame (genes by samples) 
    makeCellData <- function(expr_df) {
        genes <- rownames(expr_df)
        expr_mat = data.matrix(expr_df)
        num_cells_expressed <- (expr_mat > 0.1) + 0
        num_cells_expressed <- Matrix::rowSums(num_cells_expressed)
        fd <- data.frame(num_cells_expressed=num_cells_expressed, row.names = genes)
        fd <- new("AnnotatedDataFrame", data = fd)
        pd <- new("AnnotatedDataFrame", data = data.frame(row.names=colnames(expr_mat)))
        
        newCellDataSet(expr_mat,
            phenoData = pd,
            featureData = fd,
            lowerDetectionLimit = 0.1,
            expressionFamily = VGAM::tobit(0.1))
    }

    makeCellData3 <- function(expr_df) {
        genes <- rownames(expr_df)
        expr_mat = data.matrix(expr_df)
        num_cells_expressed <- (expr_mat > 0.1) + 0
        num_cells_expressed <- Matrix::rowSums(num_cells_expressed)
        fd <- data.frame(num_cells_expressed=num_cells_expressed, row.names = genes)
        fd <- new("AnnotatedDataFrame", data = fd)
        pd <- data.frame(row.names=colnames(expr_mat))
        # a hack to avoid error when running `partitionCells`
        pd['foo'] = 'bar'
        pd <- new("AnnotatedDataFrame", data = pd)
        
        newCellDataSet(expr_mat,
            phenoData = pd,
            featureData = fd,
            lowerDetectionLimit = 0.1)
    }

    getDEGsAsOrderingGenes <- function(cds){
        # get DEGs among clusters
        cds_expressed_genes <-  row.names(subset(fData(cds), num_cells_expressed >= 10))
        clustering_DEG_genes <- differentialGeneTest(cds[cds_expressed_genes,], 
            fullModelFormulaStr = '~Cluster',
            cores = 8)
        # order cells with top 1000 DEGs
        cds_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
        cds_ordering_genes
    }

    getHighVarGenesAsOrderingGenes <- function(cds){
        # Use genes with highest variance as ordering genes
        RowVar <- function(x, ...) {
            # from https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
            rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
        }
        # use genes with high variances for ordering cell
        gene_variances <- RowVar(exprs(cds))
        cds_ordering_genes <- names(gene_variances[order(gene_variances, decreasing = T)])[1:1000]
        cds_ordering_genes
    }

    # Run the entire Monocle-DDRTree pipeline to 
    # 1) clustering
    # 2) identify DEGs across clusters
    # 3) ordering cells/psudotime estimation
    runMonocleDDRTree <- function(cds, ordering = "de") {
        # tSNE and clustering cells
        cds <- reduceDimension(cds, 
            max_components = 2,
            norm_method = 'log',
            reduction_method = 'tSNE',
            perplexity = 5,
            verbose = T)

        n_cells <- as.numeric(dim(cds)[2])
        k <- 50 # default k for louvain clustering
        if (n_cells < 52){
            k <- n_cells - 2
        }
        cds <- clusterCells(cds, method="louvain", k=k, verbose = T)
        n_clusters <- length(unique(cds$Cluster))
        if (n_clusters > 1 && ordering == "de"){
            cds_ordering_genes <- tryCatch(
                {
                    message("Attempting to compute DEGs across clusters for ordering cells...")
                    getDEGsAsOrderingGenes(cds)
                },
                error=function(cond) {
                    message("Error encountered while computing DEGs using monocle:")
                    message(cond)
                    message("Fall back to using most variable genes for ordering cells")
                    getHighVarGenesAsOrderingGenes(cds)
                }
            )
        } else { # only 1 cluster
            message("Using most variable genes for ordering cells...")
            cds_ordering_genes <- getHighVarGenesAsOrderingGenes(cds)
        }
        cds <- setOrderingFilter(cds, ordering_genes = cds_ordering_genes)
        cds <- reduceDimension(cds, method = 'DDRTree', norm_method = 'log', 
            ncenter = NULL)
        cds <- orderCells(cds)
        return(cds)
    }


    runMonocleUMAPsimplePPT <- function(cds) {
        # 1. Noramlize and pre-process the data
        cds <- estimateSizeFactors(cds)
        cds <- preprocessCDS(cds, 
            num_dim = 50,
            norm_method = 'log',
            method = 'PCA'
            )
        # 2. Reduce the dimensionality of the data
        cds <- reduceDimension(cds, max_components = 3,
                               reduction_method = 'UMAP',
                               metric="cosine",
                               verbose = F)
        # 3. Partition the cells into supergroups
        cds <- partitionCells(cds)
        # 4. Learn the principal graph
        cds <- learnGraph(cds,
                          max_components = 3,
                          RGE_method = 'SimplePPT',
                          partition_component = T,
                          verbose = F)
        ## Not Implemented: Monocle 3 doesn't seem to support ordering by 
        # the expression of a list of genes. Will need know how to automatically
        # find the root cell instead.

        # # 5. order cells
        # cds <- orderCells(cds, 
        #     root_pr_nodes = get_correct_root_state(cds,
        #                                            cell_phenotype = 'cell_type2',
        #                                            "Multipotent progenitors"))
        return(cds)
    }

    # Convert cds object to edge_df and data_df for making plot
    # @ref: https://github.com/cole-trapnell-lab/monocle-release/blob/ea83577c511564222bd08a35a9f944b07ccd1a42/R/plotting.R#L53
    convertToDataFrames <- function(cds) {
        sample_name <- NA
        sample_state <- pData(cds)$State
        # data_dim_1 <- NA
        # data_dim_2 <- NA
        theta <- 0
        x <- 1
        y <- 2

        lib_info_with_pseudo <- pData(cds)

        reduced_dim_coords <- reducedDimK(cds)

        ica_space_df <- Matrix::t(reduced_dim_coords) %>%
          as.data.frame() %>%
          select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
          mutate(sample_name = rownames(.), sample_state = rownames(.))

        dp_mst <- minSpanningTree(cds)

        edge_df <- dp_mst %>%
          igraph::as_data_frame() %>%
          select_(source = "from", target = "to") %>%
          left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
          left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

        data_df <- t(monocle::reducedDimS(cds)) %>%
          as.data.frame() %>%
          select_(x = x, y = y) %>%
          rownames_to_column("sample_name") %>%
          mutate(sample_state) %>%
          left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

        return_rotation_mat <- function(theta) {
          theta <- theta / 180 * pi
          matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
        }
        rot_mat <- return_rotation_mat(theta)

        cn1 <- c("x", "y")
        cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
        cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
        data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
        edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
        edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
        # Drop the redundant sample_state column
        data_df[,"sample_state"] = NULL
        return(list(edge_df=edge_df, data_df=data_df))
    }


    convertToDataFrames3 <- function(cds) {
        sample_name <- NA
        sample_state <- pData(cds)$louvain_component

        x <- 1
        y <- 2
        z <- 3

        lib_info_with_pseudo <- pData(cds)

        reduced_dim_coords <- reducedDimK(cds)

        ica_space_df <- Matrix::t(reduced_dim_coords) %>%
          as.data.frame() %>%
          select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y, prin_graph_dim_3 = z) %>%
          mutate(sample_name = rownames(.), sample_state = rownames(.))

        dp_mst <- minSpanningTree(cds)

        edge_df <- dp_mst %>%
          igraph::as_data_frame() %>%
          select_(source = "from", target = "to") %>%
          left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", 
              source_prin_graph_dim_2="prin_graph_dim_2", 
              source_prin_graph_dim_3="prin_graph_dim_3"), by = "source") %>%
          left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", 
              target_prin_graph_dim_2="prin_graph_dim_2", 
              target_prin_graph_dim_3="prin_graph_dim_3"), by = "target")

        data_df <- t(monocle::reducedDimS(cds)) %>%
          as.data.frame() %>%
          select_(x = x, y = y, z = z) %>%
          rownames_to_column("sample_name") %>%
          mutate(sample_state) %>%
          left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

        # Drop the redundant sample_state column
        data_df[,"sample_state"] = NULL
        return(list(edge_df=edge_df, data_df=data_df))
    }

    #  Run the entire Monocle pipeline
    runMonoclePipeline <- function(expr_df, ordering = "de") {
        cds <- makeCellData(expr_df)
        cds <- runMonocleDDRTree(cds, ordering = ordering)
        convertToDataFrames(cds)    
    }

    runMonocle3Pipeline <- function(expr_df){
        cds <- makeCellData3(expr_df)
        cds <- runMonocleUMAPsimplePPT(cds)
        convertToDataFrames3(cds)
    }
    ''')
    runMonoclePipeline = robjects.globalenv['runMonoclePipeline']
    # Compute CPM
    rawdata = dataset.raw.to_adata().to_df().T
    data = CPM(rawdata)
    # Run Monocle
    results_monocle = runMonoclePipeline(pandas2ri.conversion.py2rpy(data), ordering=ordering)
    monocle_results = {}
    for key_idx in range(len(list(results_monocle.names))):
        key = list(results_monocle.names)[key_idx]
        df = pandas2ri.conversion.rpy2py(results_monocle[key_idx])
        monocle_results[key] = df

    monocle_results['data_df'].set_index('sample_name', inplace=True)
    monocle_results['sample_metadata'] = dataset.obs[["leiden"]].merge(
        monocle_results['data_df'],
        left_index=True,
        right_index=True
        )
    # Return
    monocle_results.update(
        {'color_by': color_by, 'plot_type': plot_type}
        )
    return monocle_results



def plot_monocle(monocle_results, debug=False):
    # Get results
    sample_metadata = monocle_results['sample_metadata']
    color_by = monocle_results.get('color_by')

    color_type = 'continuous'
    if color_by == 'State':
        color_type = 'categorical'
        
    color_column = monocle_results['sample_metadata'][color_by] if color_by else None
    sample_titles = ['<b>{}</b><br>'.format(index)+'<br>'.join('<i>{key}</i>: {value}'.format(**locals()) for key, value in rowData.items()) for index, rowData in sample_metadata.iterrows()]

    # Make a trace for the trajectory
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=1,color='#888'),
        hoverinfo='none',
        name='trajectory',
        mode='lines')

    for _, row in monocle_results['edge_df'].iterrows():
        x0, y0 = row['source_prin_graph_dim_1'], row['source_prin_graph_dim_2']
        x1, y1 = row['target_prin_graph_dim_1'], row['target_prin_graph_dim_2']
        edge_trace['x'] += (x0, x1, None)
        edge_trace['y'] += (y0, y1, None)
    
    if color_by and color_type == 'continuous':
        marker = dict(size=5, color=color_column, colorscale='Viridis', showscale=True)
        trace = go.Scatter(x=monocle_results['data_df']['x'],
            y=monocle_results['data_df']['y'],
            mode='markers',
            hoverinfo='text',
            text=sample_titles,
            marker=marker,
            name='Cells'
            )
        data = [trace, edge_trace]

    elif color_by and color_type == 'categorical' and len(color_column.unique()) <= len(s.colors):

        # Get unique categories
        unique_categories = color_column.unique()

        # Define empty list
        data = []
            
        # Loop through the unique categories
        for i, category in enumerate(unique_categories):

            # Get the color corresponding to the category
            category_color = s.colors[i]

            # Get the indices of the samples corresponding to the category
            category_indices = [i for i, sample_category in enumerate(color_column) if sample_category == category]
            
            # Create new trace
            trace = go.Scatter(x=monocle_results['data_df']['x'].values[category_indices],
                                 y=monocle_results['data_df']['y'].values[category_indices],
                                 mode='markers',
                                 hoverinfo='text',
                                 text=[sample_titles[x] for x in category_indices],
                                 name = category,
                                 marker=dict(size=5, color=category_color))
            
            # Append trace to data list
            data.append(trace)
        data.append(edge_trace)
    else:
        marker = dict(size=5)
        trace = go.Scatter(x=monocle_results['data_df']['x'],
                    y=monocle_results['data_df']['y'],
                    mode='markers',
                    hoverinfo='text',
                    text=sample_titles,
                    marker=marker)
        data = [trace, edge_trace]
    
    colored = '' if str(color_by) == 'None' else 'Colored by {}'.format(color_by)
    layout = go.Layout(title='<b>Monocle Analysis | Cell Trajectory Plot</b><br><i>{}</i>'.format(colored), hovermode='closest', margin=go.Margin(l=0,r=0,b=0,t=50), width=900,
        scene=dict(xaxis=dict(title='Component 1'), yaxis=dict(title='Component 2')))
    fig = go.Figure(data=data, layout=layout)
    fig.show()
    plt.savefig("monocle.pdf")

def run_tempora(dataset, timepoint_labels_column_name, timepoint_labels):
    robjects.r('''
    library(Seurat)
    library(Tempora)
    library(snow)
    library(RColorBrewer)
    library(igraph)

    #overwrite a function of Tempora
    BuildTrajectory <- function(object, n_pcs, difference_threshold=0.01){

      if (class(object)[1] != "Tempora"){
        stop("Not a valid Tempora object")
      }


      if (n_pcs > ncol(object@cluster.pathways.dr$rotation)){
        stop("Number of PCs selected exceeds number of PCs calculated")
      }

      significant_pathways_list <- gsva_pca <- list()
      for (i in 1:n_pcs){
        genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
        significant_pathways_list[[i]] <- object@cluster.pathways[which(rownames(object@cluster.pathways) %in% names(which(genes_scaled[,1] > 1.0 | genes_scaled[,1] < -1.0))), ]
        gsva_pca[[i]] <- colMeans(significant_pathways_list[[i]])
      }

      gsva_pca <- Reduce(rbind, gsva_pca)
      rownames(gsva_pca) <- paste0("PC", seq(1:nrow(gsva_pca)))

      mi_network <- bnlearn::aracne(as.data.frame(gsva_pca))
      edges_df <- as.data.frame(mi_network$arcs)
      edges_df$to <- as.numeric(as.character(edges_df$to))
      edges_df$from <- as.numeric(as.character(edges_df$from))
      edges_df$from_clusterscore <- unlist(sapply(edges_df$from, function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == x]))
      edges_df$to_clusterscore <- unlist(sapply(edges_df$to, function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == x]))


      edges_df$direction <- ifelse((abs(edges_df$to_clusterscore - edges_df$from_clusterscore)/(0.5*(edges_df$to_clusterscore + edges_df$from_clusterscore))) < difference_threshold, "bidirectional", "unidirectional")
      edges_df <- edges_df[-which(edges_df$from_clusterscore > edges_df$to_clusterscore), ]
      edges_df$id <- ifelse(as.numeric(edges_df$from) > as.numeric(edges_df$to), paste0(edges_df$from, edges_df$to), paste0(edges_df$to, edges_df$from))
      edges_df <- edges_df[!duplicated(edges_df$id), ]
      edges_df <- edges_df[, -6]
      edges_df$type <-  ifelse(edges_df$direction == "bidirectional", 3, 1)

      object@trajectory <- edges_df
      object@n.pcs <- n_pcs
      return(object)
    }

    #overwrite a function of Tempora
    PlotTrajectory <- function(object, plotname, layout=NULL, ...){

      if (class(object)[1] != "Tempora"){
        stop("Not a valid Tempora object")
      }

      if (is.null(object@trajectory)){
        stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
      }
      edge_graph <- igraph::graph_from_data_frame(d=object@trajectory, vertices = object@cluster.metadata, directed = T)
      if (is.null(layout)){
        l <- igraph::layout_with_sugiyama(edge_graph, layers = object@cluster.metadata$Cluster_time_score, maxiter = 1000)
        #l$layout[,2] <- 3-(rescale(object@cluster.metadata$Cluster_time_score, to=c(0,3)))
        if (length(levels(object@meta.data$Timepoints)) > 9){
          colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
          plot.igraph(edge_graph, ylim=c(-1,1), layout = l$layout, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                      vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.5, edge.width = 1.5, vertex.label.family="Arial",
                      vertex.label.color="black", edge.lty = E(edge_graph)$type, ...)
          axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
          legend(1,1, legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
        } else {
          colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
          plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = l$layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                      vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                      vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type,...)
          legend(1,1, legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
          axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
        }
        object@layouts <- l$layout

      } else {
        if (length(levels(object@meta.data$Timepoints)) > 9){
          colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
          
          plot.igraph(edge_graph, ylim=c(-1,1), layout = layout, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                      vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.5, edge.width = 1.5, vertex.label.family="Arial",
                      vertex.label.color="black", edge.lty = E(edge_graph)$type, ...)
          axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
          legend(1,1, legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
          
        } else {
          colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
          
          plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                      vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                      vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type,...)
          legend(1,1, legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
          axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
          
        }
      }
      # save
      jpeg(plotname, width=600, height=600)

      plot(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                      vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                      vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type,...)
      legend(1,1, legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      dev.off()
      
      validObject(object)
      return(object)
    }

    loadData <- function(rnaseq_data_filename, meta_data_filename, cluster_column_name, timepoint_column_name, timpoint_order_list, plotname){
        raw_counts = read.csv(rnaseq_data_filename, row.names = 1)
        meta = read.csv(meta_data_filename, row.names = 1)
        
        #create Seurat object
        # seu <- CreateSeuratObject(raw.data=raw_counts,meta.data=meta,names.delim="?",names.field=1,normalization.method="LogNormalize")
        seu <- CreateSeuratObject(counts=raw_counts,meta.data=meta,names.delim="?",names.field=1,normalization.method="LogNormalize")
        
        #convert Seurat object to Tempora object
        seu_tempora <- ImportSeuratObject(seu, assayType="RNA", clusters = cluster_column_name,
                                         timepoints = timepoint_column_name,
                                         timepoint_order = timpoint_order_list)
                                         
        seu_tempora <- CalculatePWProfiles(seu_tempora, 
                    gmt_path = "Human_GOBP_AllPathways_no_GO_iea_June_01_2020_symbol.gmt",
                    method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)
                    
                    
        #Build trajectory with 2 PCs 
        seu_tempora <- BuildTrajectory(seu_tempora, n_pcs = 2, difference_threshold = 0.01)

        #Visualize the trajectory
        seu_tempora <- PlotTrajectory(seu_tempora, plotname)
        return (seu_tempora)
    }
    ''')
    #download Human_GOBP_AllPathways_no_GO_iea_June_01_2020_symbol.gmt 
    url = "http://download.baderlab.org/EM_Genesets/June_01_2020/Human/symbol/Human_GOBP_AllPathways_no_GO_iea_June_01_2020_symbol.gmt"
    r = requests.get(url, allow_redirects=True)
    with open('Human_GOBP_AllPathways_no_GO_iea_June_01_2020_symbol.gmt', "wb") as f:
        f.write(r.content)
        
    tmp_meta_df = dataset.obs.copy()
    tmp_meta_df["Timepoints"] = tmp_meta_df[timepoint_labels_column_name]

    tmp_meta_df["leiden"] = tmp_meta_df["leiden"].astype('int')
    tmp_meta_df["leiden"] = tmp_meta_df["leiden"] + 1
    if len(tmp_meta_df["Timepoints"].unique()) > 1:
        #save preprocessed data
        tmp_meta_df.to_csv("metadata_with_clusters.csv")
        dataset.raw.to_adata().to_df().T.to_csv("expressiondata_after_preprocessing.csv")

        loadData = robjects.r['loadData']
        tempora_result = loadData("expressiondata_after_preprocessing.csv", "metadata_with_clusters.csv", "leiden", "Timepoints", timepoint_labels.split("\n"), "Tempora_plot.jpg")
