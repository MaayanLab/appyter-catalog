# %%javascript
# require.config({
#   paths: {
#     d3: 'https://cdnjs.cloudflare.com/ajax/libs/d3/5.9.2/d3',
#     jquery: 'https://code.jquery.com/jquery-3.4.1.min',
#     plotly: 'https://cdn.plot.ly/plotly-latest.min'
#   },

#   shim: {
#     plotly: {
#       deps: ['d3', 'jquery'],
#       exports: 'plotly'
#     }
#   }
# });

from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
# Basic libraries
import pandas as pd
import os
import urllib3
import requests, json
import sys
import geode
import random
from time import sleep
import time
import numpy as np
import warnings
import base64  

# Visualization
import plotly
from plotly import tools
import plotly.express as px
import plotly.graph_objs as go
plotly.offline.init_notebook_mode() # To embed plots in the output cell of the notebook

import matplotlib.pyplot as plt; plt.rcdefaults()
from matplotlib import rcParams
from matplotlib.lines import Line2D
from matplotlib_venn import venn2, venn3

import IPython
from IPython.display import HTML, display, Markdown, IFrame

import chart_studio
import chart_studio.plotly as py

# Data analysis
from itertools import combinations
import scipy.spatial.distance as dist
import scipy.stats as ss
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
from maayanlab_bioinformatics.dge.limma_voom import limma_voom_differential_expression


def printa(a=1):
  print(a)

def CPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)
        
    return data
def logCPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)
        data = np.log10(data+1)

    # Return
    return data
def log(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = data.fillna(0)
        data = np.log10(data+1)

    # Return
    return data
def qnormalization(data):
    newdata = pd.DataFrame(quantile_transform(
        data, axis=1, output_distribution='normal'))
    newdata.columns = data.columns
    newdata.index = data.index
    return newdata  

def normalize(dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization):
    normalization = 'rawdata'
    if logCPM_normalization == True:  
        data = dataset[normalization]
        normalization += '+logCPM'
        dataset[normalization] = logCPM(data)
        
    if log_normalization == True:    
        data = dataset[normalization]
        normalization += '+log'
        dataset[normalization] = log(data)
        
    if z_normalization == True:
        data = dataset[normalization]
        normalization += '+z_norm'    
        dataset[normalization] = data.T.apply(ss.zscore, axis=0).T.dropna()

    if q_normalization == True:
        data = dataset[normalization]
        normalization += '+q_norm'
        dataset[normalization] = qnormalization(data)
    return dataset, normalization

def create_download_link( df, title = "Download CSV file", filename = "data.csv"):  
    csv = df.to_csv()
    b64 = base64.b64encode(csv.encode())
    payload = b64.decode()
    html = '<a download="{filename}" href="data:text/csv;base64,{payload}" target="_blank">{title}</a>'
    html = html.format(payload=payload,title=title,filename=filename)
    return HTML(html)

def display_link(url):
    raw_html = '<a href="%s" target="_blank">%s</a>' % (url, url)
    return display(HTML(raw_html))

def run_pca(dataset, meta_id_column_name, normalization='logCPM', nr_genes=2500, color_by='auto', color_type='categorical', filter_samples=True, plot_type='interactive'):
    # Get data
    before_norm = normalization.replace("+z_norm", "").replace("+q_norm", "")
    top_genes = dataset[before_norm].var(axis=1).sort_values(ascending=False)
    
    expression_dataframe = dataset[normalization].copy()
    
    # Filter columns
    if filter_samples and dataset.get('signature_metadata'):
        selected_samples = [sample for samples in list(dataset['signature_metadata'].values())[0].values() for sample in samples]
        expression_dataframe = expression_dataframe[selected_samples]

    # Filter rows
    expression_dataframe = expression_dataframe.loc[top_genes.index[:nr_genes]]
    
    # Run PCA
    pca=PCA(n_components=3)
    pca.fit(expression_dataframe)

    # Get Variance
    var_explained = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]

    # Estimate colors
    if color_by == 'auto':

        # Add signature groups
        if dataset.get('signature_metadata'):
            A_label, B_label = list(dataset.get('signature_metadata').keys())[0].split(' vs ')
            col = []
            group_dict = list(dataset.get('signature_metadata').values())[0]
            for gsm in dataset['sample_metadata'].index:
                if gsm in group_dict['A']:
                    col.append(A_label)
                elif gsm in group_dict['B']:
                    col.append(B_label)
                else:
                    col.append('Other')
            dataset['dataset_metadata']['Sample Group'] = col
            color_by = 'Sample Group'
        else:

            # Add group column, if available
            if 'Group' in dataset['dataset_metadata'].columns:
                color_by = 'Group'
            else:
                color_by = None


    # Return
    pca_results = {'pca': pca, 'var_explained': var_explained, 
                   'dataset_metadata': dataset['dataset_metadata'][dataset['dataset_metadata'][meta_id_column_name] == expression_dataframe.columns], 
                   'color_by': color_by, 'color_type': color_type, 'nr_genes': nr_genes, 
                   'normalization': normalization, 'signature_metadata': dataset.get('signature_metadata'), 
                   'plot_type': plot_type}
    return pca_results



def plot_pca(pca_results, meta_id_column_name, meta_class_column_name,plot_type='interactive'):
    pca_transformed = pca_results['pca']
    variance_explained = pca_results['var_explained']
    meta_df = pca_results['dataset_metadata']
    
    meta_df['x'] = pca_transformed.components_[0]
    meta_df['y'] = pca_transformed.components_[1]
    meta_df['z'] = pca_transformed.components_[2]

    display(IPython.core.display.HTML('''
            <script src="/static/components/requirejs/require.js"></script>
            <script>
              requirejs.config({
                paths: {
                  base: '/static/base',
                  plotly: 'https://cdn.plot.ly/plotly-latest.min.js?noext',

                },
              });
            </script>
            '''))

    classes = meta_df[meta_class_column_name].unique().tolist()
    SYMBOLS = ['circle', 'square']
    
    if len(classes) > 10:
        def r(): return random.randint(0, 255)
        COLORS = ['#%02X%02X%02X' % (r(), r(), r())
                          for i in range(len(classes))]
    else:
        COLORS = [
            '#1f77b4',
            '#ff7f0e',
            '#2ca02c',
            '#d62728',
            '#9467bd',
            '#8c564b',
            '#e377c2',
            '#7f7f7f',
            '#bcbd22',
            '#17becf',
            ]

    data = [] # To collect all Scatter3d instances
    for (cls), meta_df_sub in meta_df.groupby([meta_class_column_name]):
        # Iteratate through samples grouped by class
        display_name = '%s' % (cls)
        # Initiate a Scatter3d instance for each group of samples specifying their coordinates
        # and displaying attributes including color, shape, size and etc.
        trace = go.Scatter3d(
            x = meta_df_sub['x'],
            y = meta_df_sub['y'],
            z = meta_df_sub['z'],

            text=meta_df_sub[meta_id_column_name],
            mode='markers',
            marker=dict(
                size=10,
                color=COLORS[classes.index(cls)], # Color by infection status
                opacity=.8,
            ),
#             name=meta_df_sub[meta_id_column_name]
            name=display_name,
        )

        data.append(trace)

    # Configs for layout and axes
    layout=dict(height=1000, width=1000, 
                title='3D PCA plot for samples',
                scene=dict(
                    xaxis=dict(title=variance_explained[0]),
                    yaxis=dict(title=variance_explained[1]),
                    zaxis=dict(title=variance_explained[2])
                    )
    )
    fig=dict(data=data, layout=layout)
    if plot_type == "interactive":
        plotly.offline.iplot(fig)
    else:
        py.image.ishow(fig)
        
        
def run_clustergrammer(dataset, meta_class_column_name, normalization='logCPM', z_score=True, nr_genes=1500, metadata_cols=None, filter_samples=True,gene_list=None):
    # Subset the expression DataFrame using top 800 genes with largest variance
    data = dataset[normalization].copy()
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
                       zip(dataset["dataset_metadata"][meta_class_column_name], expr_df_sub.columns)]
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
    display_link(clustergrammer_url)
    # Embed
    display(IPython.display.IFrame(clustergrammer_url, width="1000", height="1000"))


robjects.r('''limma <- function(rawcount_dataframe, design_dataframe, filter_genes=FALSE, adjust="BH") {
    # Load packages
    suppressMessages(require(limma))
    suppressMessages(require(edgeR))

    # Convert design matrix
    design <- as.matrix(design_dataframe)
    
    # Create DGEList object
    dge <- DGEList(counts=rawcount_dataframe)

    # Filter genes
    if (filter_genes) {
        keep <- filterByExpr(dge, design)
        dge <- dge[keep,]
    }

    # Calculate normalization factors
    dge <- calcNormFactors(dge)

    # Run VOOM
    v <- voom(dge, plot=FALSE)

    # Fit linear model
    fit <- lmFit(v, design)

    # Make contrast matrix
    cont.matrix <- makeContrasts(de=B-A, levels=design)

    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)

    # Run DE
    fit2 <- eBayes(fit2)

    # Get results
    limma_dataframe <- topTable(fit2, adjust=adjust, number=nrow(rawcount_dataframe))
    
    # Return
    results <- list("limma_dataframe"= limma_dataframe, "rownames"=rownames(limma_dataframe))
    return (results)
}
''')
def get_signatures(classes, dataset, normalization, method, meta_class_column_name, meta_id_column_name):
    expr_df = dataset['rawdata']
    meta_df = dataset["dataset_metadata"]
    signatures = dict()

    for cls1, cls2 in combinations(classes, 2):
        cls1_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls1, meta_id_column_name].tolist()
        cls2_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls2, meta_id_column_name].tolist()
        signature_label = " vs. ".join([cls1, cls2])

        if method == "limma":
            design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in cls1_sample_ids), 'B': int(x in cls2_sample_ids)} for x in expr_df.columns]).set_index('index')

            processed_data = {"expression": expr_df, 'design': design_dataframe}
            limma = robjects.r['limma']
            limma_results = pandas2ri.conversion.rpy2py(limma(pandas2ri.conversion.py2rpy(processed_data['expression']), pandas2ri.conversion.py2rpy(processed_data['design']), filter_genes=False))
            
            signature = pd.DataFrame(limma_results[0])
            signature.index = limma_results[1]
            
        elif method == "characteristic_direction":
            signature = characteristic_direction(dataset[normalization].loc[:, cls1_sample_ids], dataset[normalization].loc[:, cls2_sample_ids], calculate_sig=True)
        
        signatures[signature_label] = signature

    return signatures

# def run_cd(dataset, normalization, meta_class_column_name):
#     meta_df = dataset['dataset_metadata'].copy()
#     expr_df = dataset[normalization]
#     d_platform_cd = {} # to top up/down genes
#     cd_results = dict()

#     sample_classes = {}
#     for layout1, layout2 in combinations(meta_df[meta_class_column_name].unique(), 2):
#         sample_class = np.zeros(expr_df.shape[1], dtype=np.int32)
#         sample_class[meta_df[meta_class_column_name].values == layout1] = 1
#         sample_class[meta_df[meta_class_column_name].values == layout2] = 2
#         diff_gene_set_name = " vs. ".join([layout1, layout2])
#         sample_classes[diff_gene_set_name] = sample_class


#     for platform, sample_class in sample_classes.items():

#         cd_res = geode.chdir(expr_df.values, sample_class, expr_df.index, gamma=.5, sort=True, calculate_sig=False)
        
#         cd_coefs = np.array([x[0] for x in cd_res])
        
#         cd_coefs_series = pd.DataFrame(cd_coefs)
        
        
        
#         cd_coefs_series.index = expr_df.index
        
#         cd_coefs_series.columns = ["coef"]
#         cd_results[platform] = pd.concat([cd_coefs_series,], axis=1).sort_values("coef", ascending=False)
        
# #         cd_pvalues = np.array([x[2] for x in cd_res])
# #         cd_pvalues_series = pd.DataFrame(cd_pvalues)
# #         cd_pvalues_series.index = expr_df.index
# #         cd_pvalues_series.columns = ["significance"]
# #         cd_results[platform] = pd.concat([cd_coefs_series, cd_pvalues_series], axis=1).sort_values("coef", ascending=False)
        
        
#     return cd_results


# def get_signatures(classes, dataset, normalization, method, meta_class_column_name, meta_id_column_name):
#     expr_df = dataset['rawdata']
#     meta_df = dataset["dataset_metadata"]
#     signatures = dict()
#     if method == "limma":
#         for A, B in combinations(classes, 2):
#             group_A = meta_df.loc[meta_df[meta_class_column_name]==A, meta_id_column_name].tolist()
#             group_B = meta_df.loc[meta_df[meta_class_column_name]==B, meta_id_column_name].tolist()

#             design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in group_A), 'B': int(x in group_B)} for x in expr_df.columns]).set_index('index')

#             processed_data = {"expression": expr_df, 'design': design_dataframe}

#             limma = robjects.r['limma']
#             signature = pandas2ri.conversion.rpy2py(limma(pandas2ri.conversion.py2rpy(processed_data['expression']), pandas2ri.conversion.py2rpy(processed_data['design']), filter_genes))

#             signature_label = " vs. ".join([A, B])
#             signatures[signature_label] = signature.sort_values("t", ascending=False)
#     elif method == "characteristic_direction":
#         signatures = run_cd(dataset, normalization, meta_class_column_name)
#     return signatures


