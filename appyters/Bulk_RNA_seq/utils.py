from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
# Basic libraries
import pandas as pd
import os
import urllib3
import requests, json
import random
import time
import numpy as np
import warnings

# Visualization
import plotly
from plotly import tools
import plotly.express as px
import plotly.graph_objs as go
import matplotlib.pyplot as plt; plt.rcdefaults()
from matplotlib import rcParams

import IPython
from IPython.display import HTML, display, Markdown, IFrame


# Data analysis
from itertools import combinations
import scipy.spatial.distance as dist
import scipy.stats as ss
from sklearn.decomposition import PCA
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
import umap
from sklearn.manifold import TSNE

from plotly.offline import init_notebook_mode
init_notebook_mode(connected = False)

def check_files(fname):
    if fname == "":
        raise IOError
    if fname.endswith(".txt") == False and fname.endswith(".csv") ==False and fname.endswith(".tsv")==False:
        raise IOError
def check_df(df, col):
    if col not in df.columns:
        raise IOError
        
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
        data = np.log2(data+1)

    # Return
    return data
def log(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = data.fillna(0)
        data = np.log2(data+1)

    return data
def qnormalization(data):
  
    X_quantile_norm = quantile_normalize(data)
    return X_quantile_norm  

def normalize(dataset, current_dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization):
    normalization = current_dataset
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

def create_download_link(df, title = "Download CSV file: {}", filename = "data.csv"):  
    df.to_csv(filename)
    html = "<a href=\"./{}\" target='_blank'>{}</a>".format(filename, title.format(filename))
    return HTML(html)

def display_link(url):
    raw_html = '<a href="%s" target="_blank">%s</a>' % (url, url)
    return display(HTML(raw_html))


def display_object(counter, caption, notebook_metadata=None, saved_filename=None, df=None, istable=True, ishtml=False):
    object_info = dict()
    
    if df is not None:
        if ishtml == True:
            display(HTML(df.to_html(render_links=True, escape=False)))
        else:
            display(df)
        df.to_csv(saved_filename)
        
    if istable == True:
        if notebook_metadata is not None:
            object_info["description"] = caption
            object_info["file_name"] = saved_filename
            notebook_metadata["tables"][str(counter)] = object_info  
        

        display(Markdown("*Table {}. {}*".format(counter, caption)))
    else:
        if notebook_metadata is not None:
            object_info["description"] = caption
            object_info["file_name"] = saved_filename
            notebook_metadata["figures"][str(counter)] = object_info
        
        display(Markdown("*Figure {}. {}*".format(counter, caption)))
        
    counter += 1
    if notebook_metadata is not None:
        return counter, notebook_metadata
    else:
        return counter
    
def run_dimension_reduction(dataset, method='PCA', normalization='logCPM', nr_genes=2500, filter_samples=True, plot_type='interactive'):
    # Get data
    before_norm = normalization.replace("+z_norm", "").replace("+q_norm", "")
    top_genes = dataset[before_norm].var(axis=1).sort_values(ascending=False)
    
    expression_dataframe = dataset[normalization]
    
    # Filter columns
    if filter_samples and dataset.get('signature_metadata'):
        selected_samples = [sample for samples in list(dataset['signature_metadata'].values())[0].values() for sample in samples]
        expression_dataframe = expression_dataframe[selected_samples]

    # Filter rows
    expression_dataframe = expression_dataframe.loc[top_genes.index[:nr_genes]]
    result=None
    axis = None
    if method == 'PCA':
        # Run PCA
        pca=PCA(n_components=3)
        result = pca.fit(expression_dataframe)

        # Get Variance
        axis = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]
    elif method == "UMAP":
        u = umap.UMAP(n_components=3)
        result = u.fit_transform(expression_dataframe.T)
        axis = ['UMAP 1', 'UMAP 2', 'UMAP 3']
    elif method == "t-SNE":
        u = TSNE(n_components=3)
        result = u.fit_transform(expression_dataframe.T)
        axis = ['t-SNE 1', 't-SNE 2', 't-SNE 3']

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

    # Return
    results = {'method': method, 'result': result, 'axis': axis, 
                   'dataset_metadata': dataset['dataset_metadata'][dataset['dataset_metadata'].index == expression_dataframe.columns], 
                   'nr_genes': nr_genes, 
                   'normalization': normalization, 'signature_metadata': dataset.get('signature_metadata'), 
                   'plot_type': plot_type}
    return results


def plot_samples(pca_results, meta_class_column_name, counter, plot_name, notebook_metadata, plot_type='interactive'):
    pca_transformed = pca_results['result']
    axis = pca_results['axis']
    meta_df = pca_results['dataset_metadata']
    
    if pca_results["method"] == "PCA":
        meta_df['x'] = pca_transformed.components_[0]
        meta_df['y'] = pca_transformed.components_[1]
        meta_df['z'] = pca_transformed.components_[2]
    else:
        meta_df['x'] = [x[0] for x in pca_transformed]
        meta_df['y'] = [x[1] for x in pca_transformed]
        meta_df['z'] = [x[2] for x in pca_transformed]
        
    fig = px.scatter_3d(meta_df, x='x', y='y', z='z', color=meta_class_column_name)
    if plot_type == "interactive":
        fig.show()
    else:
        fig.show(renderer="png")
    fig.write_image(plot_name)
    
    
    caption = '3D {} plot for samples using {} genes having largest variance. The figure displays an interactive, three-dimensional scatter plot of the data. Each point represents an RNA-seq sample. Samples with similar gene expression profiles are closer in the three-dimensional space. If provided, sample groups are indicated using different colors, allowing for easier interpretation of the results.'.format(pca_results["method"], pca_results['nr_genes'])

    display(Markdown("*Figure {}. {}*".format(counter, caption)))
    object_info = dict()
    object_info["description"] = caption
    object_info["file_name"] = plot_name
    
    notebook_metadata["figures"][str(counter)] = object_info
    
    counter += 1
    return counter, notebook_metadata

def run_clustergrammer(dataset, meta_class_column_name, normalization='logCPM', z_score=True, nr_genes=1500, metadata_cols=None, filter_samples=True,gene_list=None):
    # Subset the expression DataFrame using top 800 genes with largest variance
    data = dataset[normalization]
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
    clustergrammer_url = 'https://maayanlab.cloud/clustergrammer/matrix_upload/'
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


def plot_2D_scatter(x, y, text='', title='', xlab='', ylab='', hoverinfo='text', color='black', colorscale='Blues', size=8, showscale=False, symmetric_x=False, symmetric_y=False, pad=0.5, hline=False, vline=False, return_trace=False, labels=False, plot_type='interactive', de_type='ma', plot_name="2d_plot.png"):
    range_x = [-max(abs(x))-pad, max(abs(x))+pad]if symmetric_x else None
    range_y = [-max(abs(y))-pad, max(abs(y))+pad]if symmetric_y else None
    trace = go.Scattergl(x=x, y=y, mode='markers', text=text, hoverinfo=hoverinfo, marker={'color': color, 'colorscale': colorscale, 'showscale': showscale, 'size': size})
    if return_trace:
        return trace
    else:
        if de_type == 'ma':
            annotations = [
                {'x': 1, 'y': 0.1, 'text':'<span style="color: blue; font-size: 10pt; font-weight: 600;">Down-regulated in '+labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'right', 'yanchor': 'top'},
                {'x': 1, 'y': 0.9, 'text':'<span style="color: red; font-size: 10pt; font-weight: 600;">Up-regulated in '+labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'right', 'yanchor': 'bottom'}
            ] if labels else []
        elif de_type == 'volcano':
            annotations = [
                {'x': 0.25, 'y': 1.07, 'text':'<span style="color: blue; font-size: 10pt; font-weight: 600;">Down-regulated in '+labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},
                {'x': 0.75, 'y': 1.07, 'text':'<span style="color: red; font-size: 10pt; font-weight: 600;">Up-regulated in '+labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}
            ] if labels else []
        layout = go.Layout(title=title, xaxis={'title': xlab, 'range': range_x}, yaxis={'title': ylab, 'range': range_y}, hovermode='closest', annotations=annotations)
        fig = go.Figure(data=[trace], layout=layout)
    
    if plot_type=='interactive':
        fig.show()
    else:
        fig.show(renderer="png")
    fig.write_image(plot_name)


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

def get_signatures(classes, dataset, normalization, method, meta_class_column_name, filter_genes):
    tmp_normalization = normalization.replace("+z_norm+q_norm","").replace("+z_norm","")
    raw_expr_df = dataset['rawdata']
    expr_df = dataset['rawdata']
    if filter_genes == True:
        expr_df = dataset['rawdata+filter_genes']
        
    signatures = dict()

    for cls1, cls2 in combinations(classes, 2):
        print(cls1, cls2)
        cls1_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls1, :].index.tolist() #control
        cls2_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls2,:].index.tolist() #case
        
        signature_label = " vs. ".join([cls1, cls2])
        
        if method == "limma":
            limma = robjects.r['limma']

            design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in cls1_sample_ids), 'B': int(x in cls2_sample_ids)} for x in raw_expr_df.columns]).set_index('index')

            processed_data = {"expression": raw_expr_df, 'design': design_dataframe}
            
            limma_results = pandas2ri.conversion.rpy2py(limma(pandas2ri.conversion.py2rpy(processed_data['expression']), pandas2ri.conversion.py2rpy(processed_data['design']), filter_genes=filter_genes))
                        
            signature = pd.DataFrame(limma_results[0])
            signature.index = limma_results[1]
            signature = signature.sort_values("t", ascending=False)
            
        elif method == "characteristic_direction":
            signature = characteristic_direction(dataset[tmp_normalization].loc[:, cls1_sample_ids], dataset[normalization].loc[:, cls2_sample_ids], calculate_sig=True)
            signature = signature.sort_values("CD-coefficient", ascending=False)
        elif method == "edgeR":
            edgeR = robjects.r['edgeR']
            
            edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(edgeR_results[0])
            signature.index = edgeR_results[1]
            signature = signature.sort_values("logFC", ascending=False)
        elif method == "DESeq2":
            # deseq2 receives raw counts
            DESeq2 = robjects.r['deseq2']
            DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(DESeq2_results[0])
            signature.index = DESeq2_results[1]
            signature = signature.sort_values("log2FoldChange", ascending=False)
                        
            
        signatures[signature_label] = signature

    return signatures        


def run_volcano(signature, signature_label, dataset, pvalue_threshold, logfc_threshold, plot_type):
    color = []
    text = []
    for index, rowData in signature.iterrows():
        if "AveExpr" in rowData.index: # limma
            expr_colname = "AveExpr"
            pval_colname = "P.Value"
            logfc_colname = "logFC"
        elif "logCPM" in rowData.index: #edgeR
            expr_colname = "logCPM"
            pval_colname = "PValue"
            logfc_colname = "logFC"
        elif "baseMean" in rowData.index: #DESeq2
            expr_colname = "baseMean"
            pval_colname = "pvalue"
            logfc_colname = "log2FoldChange"
        # Text
        text.append('<b>'+index+'</b><br>Avg Expression = '+str(round(rowData[expr_colname], ndigits=2))+'<br>logFC = '+str(round(rowData[logfc_colname], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData[pval_colname])+'<br>FDR = '+'{:.2e}'.format(rowData[pval_colname]))

        # Color
        if rowData[pval_colname] < pvalue_threshold:
            if rowData[logfc_colname] < -logfc_threshold:
                color.append('blue')
            elif rowData[logfc_colname] > logfc_threshold:
                color.append('red')
            else:
                color.append('black')

        else:
            color.append('black')

    volcano_plot_results = {'x': signature[logfc_colname], 'y': -np.log10(signature[pval_colname]), 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results
        

def plot_volcano(volcano_plot_results):
    spacer = ' '*50
    plot_name = "volcano_plot_{}.png".format(volcano_plot_results['signature_label'])
    plot_2D_scatter(
        x=volcano_plot_results['x'],
        y=volcano_plot_results['y'],
        text=volcano_plot_results['text'],
        color=volcano_plot_results['color'],
        symmetric_x=True,
        xlab='log2FC',
        ylab='-log10P',
        title='<b>{volcano_plot_results[signature_label]} Signature | Volcano Plot</b>'.format(**locals()),
        labels=volcano_plot_results['signature_label'].split(' vs. '),
        plot_type=volcano_plot_results['plot_type'],
        de_type='volcano',
        plot_name=plot_name
    )        
    return plot_name

def run_maplot(signature, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive'):

    # Loop through signature
    color = []
    text = []
    for index, rowData in signature.iterrows():
        if "AveExpr" in rowData.index: # limma
            expr_colname = "AveExpr"
            pval_colname = "adj.P.Val"
            logfc_colname = "logFC"
        elif "logCPM" in rowData.index: #edgeR
            expr_colname = "logCPM"
            pval_colname = "PValue"
            logfc_colname = "logFC"
        elif "baseMean" in rowData.index: #DESeq2
            expr_colname = "baseMean"
            pval_colname = "padj"
            logfc_colname = "log2FoldChange"
        # Text
        text.append('<b>'+index+'</b><br>Avg Expression = '+str(round(rowData[expr_colname], ndigits=2))+'<br>logFC = '+str(round(rowData[logfc_colname], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData[pval_colname])+'<br>FDR = '+'{:.2e}'.format(rowData[pval_colname]))

        # Color
        if rowData[pval_colname] < pvalue_threshold:
            if rowData[logfc_colname] < -logfc_threshold:
                color.append('blue')
            elif rowData[logfc_colname] > logfc_threshold:
                color.append('red')
            else:
                color.append('black')

        else:
            color.append('black')
    
    # Return 
    volcano_plot_results = {'x': signature[expr_colname], 'y': signature[logfc_colname], 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results

def plot_maplot(volcano_plot_results):    
    plot_name = "ma_plot_{}.png".format(volcano_plot_results['signature_label'])
    plot_2D_scatter(
        x=volcano_plot_results['x'],
        y=volcano_plot_results['y'],
        text=volcano_plot_results['text'],
        color=volcano_plot_results['color'],
        symmetric_y=True,
        xlab='Average Expression',
        ylab='logFC',
        title='<b>{volcano_plot_results[signature_label]} Signature | MA Plot</b>'.format(**locals()),
        labels=volcano_plot_results['signature_label'].split(' vs. '),
        plot_type=volcano_plot_results['plot_type'],
        plot_name=plot_name
        
    )
    return plot_name


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

def submit_enrichr_geneset(geneset, label=''):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
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
            'GO_Molecular_Function_'+go_version: 'Gene Ontology Molecular Function ('+go_version+' version)',
            'GO_Cellular_Component_'+go_version: 'Gene Ontology Cellular Component ('+go_version+' version)'
        }
    elif library_type == "pathway":
        # Libraries
        libraries = {
            'KEGG_2019_Human': 'KEGG Pathways',
            'WikiPathways_2019_Human': 'WikiPathways',
            'Reactome_2016': 'Reactome Pathways'
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
            'ENCODE_TF_ChIP-seq_2015': 'ENCODE (experimentally validated targets)',
            'ARCHS4_TFs_Coexp': 'ARCHS4 (coexpressed genes)'
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
    

def plot_library_barchart(enrichr_results, gene_set_library, signature_label, sort_results_by='pvalue', nr_genesets=15, height=400, plot_type='interactive', plot_name="barchart.png"):
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
    labels = signature_label.split(' vs. ')
    annotations = [
        {'x': 0.25, 'y': 1.06, 'text': '<span style="color: #FA8072; font-size: 10pt; font-weight: 600;">Up-regulated in ' +
            labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},
        {'x': 0.75, 'y': 1.06, 'text': '<span style="color: #87CEFA; font-size: 10pt; font-weight: 600;">Down-regulated in ' +
            labels[-1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}
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
    if plot_type=='interactive':
        fig.show()
    else:
        fig.show(renderer='png')
    fig.write_image(plot_name)




def results_table(enrichment_dataframe, source_label, target_label, notebook_metadata, table_counter):

    # Get libraries
    for gene_set_library in enrichment_dataframe['gene_set_library'].unique():

        # Get subset
        enrichment_dataframe_subset = enrichment_dataframe[enrichment_dataframe['gene_set_library'] == gene_set_library].copy()

        # Get unique values from source column
        enrichment_dataframe_subset[source_label] = [x.split('_')[0] for x in enrichment_dataframe_subset['term_name']]
        enrichment_dataframe_subset = enrichment_dataframe_subset.sort_values(['FDR', 'pvalue']).rename(columns={'pvalue': 'P-value'}).drop_duplicates(source_label)

        # Add links and bold for significant results
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
        if source_label == "Transcription Factor":
            additional_description = ". The table contains scrollable tables displaying the results of the Transcription Factor (TF) enrichment analysis generated using Enrichr. Every row represents a TF; significant TFs are highlighted in bold."
        elif source_label == "Kinase":
            additional_description = ". The table contains browsable tables displaying the results of the Protein Kinase (PK) enrichment analysis generated using Enrichr. Every row represents a PK; significant PKs are highlighted in bold."    
        elif source_label == "miRNA":
            additional_description = ". The figure contains browsable tables displaying the results of the miRNA enrichment analysis generated using Enrichr. Every row represents a miRNA; significant miRNAs are highlighted in bold."
        filename="Enrichment_analysis_{}_{}.csv".format(source_label, gene_set_library)
        table_counter, notebook_metadata = display_object(table_counter, gene_set_library+additional_description, notebook_metadata, saved_filename=filename, istable=True)
        display(create_download_link(enrichment_dataframe_subset, filename=filename))
        
    return table_counter, notebook_metadata

def display_table(analysis_results, source_label, notebook_metadata, table_counter):
    
    # Plot Table
    return results_table(analysis_results['enrichment_dataframe'].copy(), source_label=source_label, notebook_metadata=notebook_metadata, target_label='target', table_counter=table_counter)



def run_l1000cds2(signature, nr_genes=500, signature_label='', plot_type='interactive'):
    # Define results
    l1000cds2_results = {'signature_label': signature_label}

    # Define upperGenes Function
    upperGenes = lambda genes: [gene.upper() for gene in genes]

    # Get Data
    data = {"upGenes":upperGenes(signature.index[:nr_genes]),"dnGenes":upperGenes(signature.index[-nr_genes:])}

    # Loop through aggravate:
    for aggravate in [True, False]:

        # Send to API
        config = {"aggravate":aggravate,"searchMethod":"geneSet","share":True,"combination":False,"db-version":"latest"}
        r = requests.post('http://maayanlab.cloud/L1000CDS2/query',data=json.dumps({"data":data,"config":config}),headers={'content-type':'application/json'})
        label = 'mimic' if aggravate else 'reverse'

        # Add results
        resGeneSet = r.json()
        if resGeneSet.get('topMeta'):
            l1000cds2_dataframe = pd.DataFrame(resGeneSet['topMeta'])[['cell_id', 'pert_desc', 'pert_dose', 'pert_dose_unit', 'pert_id', 'pert_time', 'pert_time_unit', 'pubchem_id', 'score', 'sig_id']].replace('-666', np.nan)
            l1000cds2_results[label] = {'url': 'http://maayanlab.cloud/L1000CDS2/#/result/{}'.format(resGeneSet['shareId']), 'table': l1000cds2_dataframe}
        else:
            l1000cds2_results[label] = None
    l1000cds2_results['plot_type'] = plot_type
    # Return
    return l1000cds2_results
    
def plot_l1000cds2(l1000cds2_results, counter, notebook_metadata, nr_drugs=7, height=300, plot_name="L1000CDS2.png"):

    # Check if there are results
    if not l1000cds2_results['mimic'] or not l1000cds2_results['reverse']:
        print('### No results were found.\n This is likely due to the fact that the gene identifiers were not recognized by L1000CDS<sup>2</sup>. Please note that L1000CDS<sup>2</sup> currently only supports HGNC gene symbols (https://www.genenames.org/). If your dataset uses other gene identifier systems, such as Ensembl IDs or Entrez IDs, consider converting them to HGNC. Automated gene identifier conversion is currently under development.')
        
    else:
        # Bar charts
        fig = tools.make_subplots(rows=1, cols=2, print_grid=False);
        for i, direction in enumerate(['mimic', 'reverse']):
            drug_counts = l1000cds2_results[direction]['table'].groupby('pert_desc').size().sort_values(ascending=False).iloc[:nr_drugs].iloc[::-1]

            # Get Bar
            bar = go.Bar(
                x=drug_counts.values,
                y=drug_counts.index,
                orientation='h',
                name=direction.title(),
                hovertext=drug_counts.index,
                hoverinfo='text',
                marker={'color': '#FF7F50' if direction=='mimic' else '#9370DB'}
            )
            fig.append_trace(bar, 1, i+1)
            
            # Get text
            text = go.Scatter(
                x=[max(bar['x'])/50 for x in range(len(bar['y']))],
                y=bar['y'],
                mode='text',
                hoverinfo='none',
                showlegend=False,
                text=drug_counts.index,
                textposition="middle right",
                textfont={'color': 'black'}
            )
            fig.append_trace(text, 1, i+1)

        fig['layout'].update(height=height, title='<b>L1000CDS<sup>2</sup> {}| Small Molecule Query</b><br><i>Top small molecules</i>'.format(l1000cds2_results['signature_label']), hovermode='closest')
        fig['layout']['xaxis1'].update(domain=[0,0.5])
        fig['layout']['xaxis1'].update(title='<br>Count')
        fig['layout']['xaxis2'].update(title='<br>Count')
        fig['layout']['xaxis2'].update(domain=[0.5,1])
        fig['layout']['yaxis1'].update(showticklabels=False)
        fig['layout']['yaxis2'].update(showticklabels=False)
        fig['layout']['margin'].update(l=10, t=95, r=0, b=45, pad=5)

        if l1000cds2_results['plot_type'] == 'interactive':
            fig.show()       
        else:
            fig.show(renderer="png")
        
        counter, notebook_metadata = display_object(counter, "Top {} Mimic/Reverse Small Molecule from L1000CDS2 for {}.".format(nr_drugs, l1000cds2_results['signature_label']), notebook_metadata, saved_filename=plot_name, istable=False)
        notebook_metadata["figures"][str(counter-1)]["file_desc"] = [{"name": "Mimic Signature Query Results", "link":l1000cds2_results['mimic']['url']},
                                                                     {"name": "Reverse Signature Query Results", "link":l1000cds2_results['reverse']['url']},]
        
        # Links        
        display(Markdown(' *Mimic Signature Query Results*:'))
        display_link(l1000cds2_results['mimic']['url'])
        display(Markdown(' *Reverse Signature Query Results*:'))
        display_link(l1000cds2_results['reverse']['url'])
        
    return counter, notebook_metadata
        
def run_l1000fwd(signature, nr_genes=500, signature_label=''):
    # Define results
    l1000fwd_results = {'signature_label': signature_label}

    # Define upperGenes Function
    upperGenes = lambda genes: [gene.upper() for gene in genes]

    # Get Data
    payload = {"up_genes":upperGenes(signature.index[:nr_genes]),"down_genes":upperGenes(signature.index[-nr_genes:])}

    # Get URL
    L1000FWD_URL = 'https://maayanlab.cloud/L1000FWD/'

    # Get result
    response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
    if 'KeyError' in response.text:
        l1000fwd_results['result_url'] = None
    else:
        # Get ID and URL
        result_id = response.json()['result_id']
        l1000fwd_results['result_url'] = 'https://maayanlab.cloud/l1000fwd/vanilla/result/'+result_id
        l1000fwd_results['result_id'] = result_id

        # Get Top
        l1000fwd_results['signatures'] = requests.get(L1000FWD_URL + 'result/topn/' + result_id).json()


    # Return
    return l1000fwd_results
    
def plot_l1000fwd(l1000fwd_results, figure_counter, table_counter, notebook_metadata, nr_drugs=7, height=300):
    
    # Check if results
    if l1000fwd_results['result_url']:

        # Display IFrame
        display_link(l1000fwd_results['result_url'])
        display(IFrame(l1000fwd_results['result_url'], width="1000", height="1000"))
        figure_counter, notebook_metadata = display_object(figure_counter, "L1000FWD plot for visualizing top mimickers (red) and top reversers (blue) of {}.".format(l1000fwd_results["signature_label"]), notebook_metadata, saved_filename=l1000fwd_results['result_url'], istable=False)
        # Display tables
        for direction, signature_list in l1000fwd_results['signatures'].items():

            # Fix dataframe
            rename_dict = {'sig_id': 'Signature ID', 'pvals': 'P-value', 'qvals': 'FDR', 'zscores': 'Z-score', 'combined_scores': 'Combined Score'}
            signature_dataframe = pd.DataFrame(signature_list)[list(rename_dict.keys())].rename(columns=rename_dict).sort_values('P-value').rename_axis('Rank')
            signature_dataframe.index = [x + 1 for x in range(len(signature_dataframe.index))]
            signature_txt = signature_dataframe.to_csv(sep='\t')

            # Display table
            pd.set_option('max.colwidth', -1)
            signature_dataframe['Signature ID'] = ['<a href="http://maayanlab.cloud/dmoa/sig/{x}" target="_blank">{x}</a>'.format(**locals()) for x in signature_dataframe['Signature ID']]
            table_html = signature_dataframe.to_html(escape=False, classes='w-100')
            
            display(HTML('<style>.w-100{{width: 100% !important;}}</style><div style="max-height: 250px; overflow-y: auto; margin-bottom: 25px;">{}</div>'.format(table_html)))
            filename = "{} Signatures for {}.csv".format(direction.title(), l1000fwd_results["signature_label"])
            table_counter, notebook_metadata = display_object(table_counter, "L1000FWD Results: {} Signatures of {}".format(direction.title(), l1000fwd_results["signature_label"]), notebook_metadata, saved_filename=filename, istable=True)
            display(create_download_link(signature_dataframe, filename=filename))
            
    return figure_counter, table_counter, notebook_metadata

