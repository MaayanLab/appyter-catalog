# Basic libraries
import pandas as pd
import os
import urllib3
import requests, json
import sys
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

import IPython
from IPython.display import HTML, display, Markdown, IFrame

import chart_studio
import chart_studio.plotly as py


def create_download_link(df, title = "Download CSV file: {}", filename = "data.csv"):  
    df.to_csv(filename)
    html = "<a href=\"./{}\" target='_blank'>{}</a>".format(filename, title.format(filename))
    return HTML(html)

def display_link(url):
    raw_html = '<a href="%s" target="_blank">%s</a>' % (url, url)
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


def submit_enrichr_geneset(geneset, label=''):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
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

def run_enrichr(geneset, signature_label):
    # Submit to Enrichr
    enrichr_ids = {"result": submit_enrichr_geneset(geneset=geneset, label=signature_label)}
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
            'GO_Molecular_Function_'+go_version: 'Gene Ontology Molecular Function ('+go_version+' version)',
            'GO_Cellular_Component_'+go_version: 'Gene Ontology Cellular Component ('+go_version+' version)'
        }
    elif library_type == "pathway":
        # Libraries
        libraries = {
            'KEGG_2016': 'KEGG Pathways',
            'WikiPathways_2016': 'WikiPathways',
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
            'ChEA_2016': 'A. ChEA (experimentally validated targets)',
            'ENCODE_TF_ChIP-seq_2015': 'B. ENCODE (experimentally validated targets)',
            'ARCHS4_TFs_Coexp': 'C. ARCHS4 (coexpressed genes)'
        }
    elif library_type == "ke":
        # Libraries
        libraries = {
            'KEA_2015': 'A. KEA (experimentally validated targets)',
            'ARCHS4_Kinases_Coexp': 'B. ARCHS4 (coexpressed genes)'
        }
    elif library_type == "mirna":
        libraries = {
        'TargetScan_microRNA_2017': 'A. TargetScan (experimentally validated targets)',
        'miRTarBase_2017': 'B. miRTarBase (experimentally validated targets)'
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
    

def plot_library_barchart(enrichr_results, gene_set_library, signature_label, sort_results_by='pvalue', nr_genesets=15, height=400, plot_type='interactive'):
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
        plotly.offline.iplot(fig)
    else:
        py.image.ishow(fig)



def display_result_table(enrichment_dataframe, table_name, table_counter):
    source_label = "TF"
    # Add links
    enrichment_dataframe[source_label] = ['<a href="http://amp.pharm.mssm.edu/Harmonizome/gene/{}" target="_blank">{}</a>'.format(x.split(" ")[0], x)for x in enrichment_dataframe[source_label]]
   

    # Add rank
    enrichment_dataframe['Rank'] = ['<b>'+str(x+1)+'</b>' for x in range(len(enrichment_dataframe.index))]
    
    # Bold for significant results
    if 'FDR' in enrichment_dataframe.columns:
        enrichment_dataframe[source_label] = [rowData[source_label].replace('target="_blank">', 'target="_blank"><b>').replace('</a>', '*</b></a>') if float(rowData['FDR']) < 0.05 else rowData[source_label] for index, rowData in enrichment_dataframe.iterrows()]
    
    
    # Add overlapping genes with tooltip
    enrichment_dataframe['nr_overlapping_genes'] = [len(x) for x in enrichment_dataframe['Overlapping_Genes']]

    
    # Convert to HTML
    pd.set_option('max.colwidth', -1)
    if 'FDR' in enrichment_dataframe.columns:
        html_table = enrichment_dataframe.head(50)[['Rank', source_label, 'FDR', 'FET p-value', 'Library']].to_html(escape=False, index=False, classes='w-100')
    else:
        html_table = enrichment_dataframe.head(50)[['Rank', source_label, 'Score', 'Library']].to_html(escape=False, index=False, classes='w-100')
    html_results = '<div style="max-height: 200px; overflow-y: scroll;">{}</div>'.format(html_table)
    
    table_counter += 1
    # Display table
    display(HTML(html_results))
    
    additional_description = "{} library result from ChEA3. The table contains scrollable tables displaying the results of the Transcription Factor (TF) enrichment analysis from ChEA3. Every row represents a TF; significant TFs are highlighted in bold.".format(table_name)
    
    
    display_object(table_counter, additional_description, istable=True)
    display(create_download_link(enrichment_dataframe, filename="ChEA3_Enrichment_analysis_{}.csv".format(table_name)))
    return table_counter

def display_table(analysis_results, source_label, table_counter):
    
    # Plot Table
    return results_table(analysis_results['enrichment_dataframe'].copy(), source_label=source_label, target_label='target', table_counter=table_counter)


def run_chea3(signature, signature_label=''):
    chea3_URL = "https://maayanlab.cloud/chea3/api/enrich/"
    chea3_query = {"query_name": signature_label, "gene_set": signature}
    # Get result
    response = requests.post(chea3_URL, json=chea3_query)
    if 'KeyError' in response.text:
        chea3_result = None
    else:
        # Get ID and URL
        return response.json()
      
