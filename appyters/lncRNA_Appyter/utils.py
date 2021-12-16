
import requests
import os,json
import urllib.request
import ssl
import re
import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from plotly.offline import iplot
import plotly.graph_objs as go
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, CustomJS, ColumnDataSource, Select, Legend, Paragraph, LinearColorMapper, ColorBar, CategoricalColorMapper
from bokeh.layouts import row, column
from bokeh.palettes import Pastel2, Set2, Set1, Colorblind
from bokeh.models import Arrow, NormalHead
import hashlib
import colorcet as cc
import itertools
import random
from pyvis.network import Network
from scipy.sparse import load_npz
import copy
from bokeh.palettes import Colorblind, Bokeh, Category20,Category20b, Accent
import s3fs
import seaborn as sns
import matplotlib.pyplot as plt
from bokeh.models import Title

# Get Enrichr link
def Enrichr_API(enrichr_gene_list, description):

    short_id = ''

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(enrichr_gene_list)
    description = description
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    short_id = data["shortId"]
    return('https://maayanlab.cloud/Enrichr/enrich?dataset='+ str(short_id))

# Load Enrichr libraries
def loadLibrary(library: str, overwrite: bool = False) -> str:
    ssl._create_default_https_context = ssl._create_unverified_context
    if not os.path.exists("gmts/"+library +'.gmt' or overwrite):
        os.makedirs("gmts", exist_ok=True)
        print("Download Enrichr geneset library")
        urllib.request.urlretrieve("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="+library, "gmts/"+library+".gmt")
    else:
        print("File cached. To reload use loadLibrary(\""+library+"\", overwrite=True) instead.")
    return("gmts/"+library+".gmt")

# Predict functions based on mean pearson correaltion for each term in a library 
def predict_functions(library, matrix, query):
    library_path = loadLibrary(library)
    open_gmt = open(library_path,'r')
    library_dict = {}
    for line in open_gmt.readlines():
        line = line.strip().split('\t')
        term = line[0]
        gene_set = line[2:]
        library_dict[term]=gene_set
    open_gmt.close()  

    all_terms = []
    all_scores = []

    for lib_term, gene_set in library_dict.items():
        all_terms.append(lib_term)
        lib_term_set = list(set(gene_set)&set(matrix.index))
        lib_term_set = [x for x in lib_term_set if x!= query]
        all_scores.append(np.mean(matrix.loc[lib_term_set]["Pearson's Correlation Coefficient"]))

    df_results = pd.DataFrame({'Term':all_terms,'Mean Pearson Correlation':all_scores})
    df_results = df_results.sort_values(by ='Mean Pearson Correlation',ascending=False)
    return(df_results)

def plot_bar(df,title,x_label,y_label,filename):
    df_dict = {'x':list(df.index),'y':df[list(df.columns)[0]]}
    bar = figure(x_range=df_dict['x'], height=500, width=1000, title=title,toolbar_location="right", tools=["hover",'save'],tooltips=[ (x_label, "@x"),(y_label, "@y")])
    bar.vbar(x='x', top='y', width=0.7, source=df_dict)
    bar.xgrid.grid_line_color = None
    #bar.y_range.start = 0
    bar.xaxis.major_label_orientation = 1
    bar.xaxis.axis_label = x_label
    bar.yaxis.axis_label = y_label
    show(bar)

    # Save image
    plt.subplots(figsize=(20, 8))
    sns.barplot(x=df_dict['x'], y=df_dict['y'],color='steelblue')
    plt.ylabel(y_label,fontsize=20)
    plt.xlabel(x_label,fontsize=20)
    plt.xticks(rotation=60,fontsize=15,ha='right')
    plt.title(title,fontsize=20)
    plt.savefig(filename+'.png',facecolor='white',bbox_inches='tight')
    plt.close()
  

# Plot the top terms for each prediction library
def plot_results(library_names, results_dfs, file_name, top_results=15):
    
    fig = make_subplots(rows=1, cols=2, print_grid=False,shared_xaxes=False)
    max_scores = []
    for i in range(0,2):
        results_df = results_dfs[i][0:top_results].sort_values(by='Mean Pearson Correlation')
        library_name = library_names[i]
        max_scores.append(np.max(results_df['Mean Pearson Correlation']))
        bar = go.Bar(x=results_df['Mean Pearson Correlation'],
            y=results_df['Term'],
            orientation='h',
            name=library_name,
            showlegend=False,
            hovertext=['<b>Term: {Term}</b><br><b>Mean Pearson Correlation</b>: <i>{Mean Pearson Correlation:.3}</i>'.format(**rowData) for index, rowData in results_df[0:top_results].iterrows()],
            hoverinfo='text', 
            marker={'color': 'lightskyblue'})
        fig.append_trace(bar, 1, i+1)
        
        #Get text
        text_shortened = ['<b>{}</b>'.format(rowData['Term']) for index, rowData in results_df[0:top_results].iterrows()]
        text_shortened = [str(x[0:50] + '...' + x[-25::]) if len(x) > 75 else x for x in text_shortened] # if text is longer than image, shorten text
        text = go.Scatter(
            x=[max(bar['x'])/50 for x in range(len(bar['y']))],
            y=bar['y'],
            mode='text',
            hoverinfo='none',
            showlegend=False,
            text=text_shortened,
            textposition="middle right",
            textfont={'color': 'black','size':8})
        fig.append_trace(text, 1, i+1)
    
    annotations= [{'x': 0.25, 'y': 1.1, 'text': '<span style="color: black; font-size: 12pt; font-weight: 600;">' +library_names[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},{'x': 0.75, 'y': 1.1, 'text': '<span style="color: black; font-size: 12pt; font-weight: 600;">' +library_names[1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}]
    fig['layout'].update(height = 500, hovermode='closest', annotations=annotations)
    fig.update_layout(title='',height = 500,title_font_size = 25,title_x=0.5)
    
    fig['layout']['xaxis1'].update(domain=[0, 0.49], title='Mean Pearson Correlation' ,range=(0,max_scores[0]+max_scores[0]*.01))
    fig['layout']['xaxis2'].update(domain=[0.51, 1], title='Mean Pearson Correlation',range=(0,max_scores[1]+max_scores[1]*.01))
    fig['layout']['yaxis1'].update(showticklabels=False)
    fig['layout']['yaxis2'].update(showticklabels=False)
    fig['layout']['margin'].update(l=30, t=65, r=30, b=35)

    # Show plot
    iplot(fig)

    fig = make_subplots(rows=1, cols=2, print_grid=False,shared_xaxes=False)
    max_scores = []
    for i in range(0,2):
        results_df = results_dfs[i][0:top_results].sort_values(by='Mean Pearson Correlation')
        library_name = library_names[i]
        max_scores.append(np.max(results_df['Mean Pearson Correlation']))
        bar = go.Bar(x=results_df['Mean Pearson Correlation'],
            y=results_df['Term'],
            orientation='h',
            name=library_name,
            showlegend=False,
            hovertext=['<b>Term: {Term}</b><br><b>Mean Pearson Correlation</b>: <i>{Mean Pearson Correlation:.3}</i>'.format(**rowData) for index, rowData in results_df[0:top_results].iterrows()],
            hoverinfo='text', 
            marker={'color': 'lightskyblue'})
        fig.append_trace(bar, 1, i+1)
        
        #Get text
        text_shortened = ['<b>{}</b>'.format(rowData['Term']) for index, rowData in results_df[0:top_results].iterrows()]
        text_shortened = [str(x[0:50] + '...' + x[-25::]) if len(x) > 75 else x for x in text_shortened] # if text is longer than image, shorten text
        text = go.Scatter(
            x=[max(bar['x'])/50 for x in range(len(bar['y']))],
            y=bar['y'],
            mode='text',
            hoverinfo='none',
            showlegend=False,
            text=text_shortened,
            textposition="middle right",
            textfont={'color': 'black','size':10})
        fig.append_trace(text, 1, i+1)
    
    annotations= [{'x': 0.25, 'y': 1.05, 'text': '<span style="color: black; font-size: 12pt; font-weight: 600;">' +library_names[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},{'x': 0.75, 'y': 1.05, 'text': '<span style="color: black; font-size: 12pt; font-weight: 600;">' +library_names[1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}]
    fig['layout'].update(height = 600, width = 1000, hovermode='closest', annotations=annotations)
    fig.update_layout(title='',title_font_size = 30,title_x=0.5)
    
    fig['layout']['xaxis1'].update(domain=[0, 0.49], title='Mean Pearson Correlation' ,range=(0,max_scores[0]+max_scores[0]*.01))
    fig['layout']['xaxis2'].update(domain=[0.51, 1], title='Mean Pearson Correlation',range=(0,max_scores[1]+max_scores[1]*.01))
    fig['layout']['yaxis1'].update(showticklabels=False)
    fig['layout']['yaxis2'].update(showticklabels=False)
    fig['layout']['margin'].update(l=30, t=65, r=30, b=35)    
    fig.write_image(file_name+".png")
    #fig.write_image(file_name+".svg")

def str_to_int(string, mod):
    string = re.sub(r"\([^()]*\)", "", string).strip()
    byte_string = bytearray(string, "utf8")
    return int(hashlib.sha256(byte_string).hexdigest(), base=16)%mod

def plot_scatter(x,y,values,query,title,min_val,max_val,arrow_loc,rank,filename):
    fig = plt.figure(figsize=(20,15))
    plt.scatter(x, y, s=6, c=values, cmap=plt.cm.get_cmap('seismic'), vmin = min_val, vmax = max_val)
    cb = plt.colorbar()
    cb.set_label(label='Z-score',fontsize=15)
    plt.xlabel('UMAP 1',fontsize=20)
    plt.ylabel('UMAP 2',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.annotate('', (arrow_loc['x'], arrow_loc['y']),xytext=(arrow_loc['x']+0.2, arrow_loc['y']+0.6),arrowprops=dict(facecolor='black', headwidth= 8, width=0.1))
    plt.title(title,fontsize=25)
    plt.savefig(filename+query+'_'+title+'_rank'+str(rank) +'.png', bbox_inches='tight', dpi=300, facecolor='white')
    #plt.savefig(filename+query+'_'+title+'_rank'+str(rank) +'.svg', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()


def plot_dynamic_scatter(umap_df, values_dict, option_list, sample_names, caption_text, category_list_dict=None, location='right', category=True, dropdown=False, figure_counter=0, additional_info=None, color_by_title='', first_selection = None, highlight_query=None, static_images_save = [], file_path=''):
    if first_selection == None:
        first_selection = option_list[0]

    #Find the min and max scores across all options
    
    max_val = np.max([np.max(values_dict[option]) for option in values_dict.keys()])
    min_val = np.min([np.min(values_dict[option]) for option in values_dict.keys()])
    max_val_round = np.ceil(np.max([np.abs(max_val),np.abs(min_val)]))
    max_val = max_val_round
    min_val = np.negative(max_val_round)

    # init plot 
    source = ColumnDataSource(data=dict(x=umap_df["x"], y=umap_df["y"], values=values_dict[first_selection], names=sample_names))
    
    # node size
    if umap_df.shape[0] > 1000:
        node_size = 4
    else:
        node_size = 6
    plot = figure(plot_width=1000, plot_height=800,sizing_mode='scale_both')      
 
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
            if key in ['all tissues','all cell lines']:
                palette = list(itertools.chain(*zip(Pastel2[8], Set2[8], Set1[9], Colorblind[8])))
                unused_color = list(itertools.chain(*zip(Pastel2[8], Set2[8], Set1[9], Colorblind[8])))
                factors_dict[key] = category_list_dict[key]
                colors_dict[key] = list()
                for category_name in factors_dict[key]:
                    color_for_category = palette[str_to_int(category_name, len(palette))]
                    
                    if color_for_category not in unused_color:
                        if len(unused_color) > 0:
                            color_for_category = unused_color[0]                        
                        else:
                            color_for_category = random.sample(palette,1)[0]
                    
                    colors_dict[key].append(color_for_category)
                    if color_for_category in unused_color:
                        unused_color.remove(color_for_category)
            else:
                factors_dict[key] = category_list_dict[key]
                colors_dict[key] = list()
                for category_name in factors_dict[key]:
                    if category_name == 'other':
                        colors_dict[key].append('#cccccc')
                    else:
                        colors_dict[key].append('#e41a1c')
           
        color_mapper = CategoricalColorMapper(factors=factors_dict[first_selection], palette=colors_dict[first_selection])
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
        unique_category_dict = dict()
        for option in option_list:
            unique_category_dict[option] = sorted(list(set(values_dict[option])))
        
        # map category to color
        # color is mapped by its category name 
        # if a color is used by other categories, use another color
        colors_dict = dict()
        for key in values_dict.keys():
            colors_dict[key]= cc.CET_D1A

        color_mapper = LinearColorMapper(palette=colors_dict[first_selection] , low=min_val, high=max_val)
        #color_mapper = LinearColorMapper(palette=colors_dict[first_selection] , low=min(values_dict[first_selection]), high=max(values_dict[first_selection]))
        color_bar = ColorBar(color_mapper=color_mapper)
        plot.add_layout(color_bar, 'right')
        scatter = plot.scatter('x', 'y', size=node_size,  source=source, color={'field': 'values', 'transform': color_mapper})
       
    if additional_info is not None:
            tooltips = [
            ("lncRNA", "@names"),
            ("Label", "@values"),
            ("p-value", "@info")
        ]
    else:
        tooltips = [
            ("lncRNA", "@names"),
            ("Label", "@values"),
        ]
    
    if highlight_query!=None:
        arrow_loc = umap_df.loc[highlight_query]
        plot.add_layout(Arrow(end=NormalHead(size=10,fill_color="black"),x_start=arrow_loc['x']+0.2, y_start=arrow_loc['y']+0.6, x_end=arrow_loc['x'], y_end=arrow_loc['y']))
    
    plot.add_tools(HoverTool(tooltips=tooltips))
    plot.output_backend = "webgl"
    
    plot.xaxis.axis_label = "UMAP 1"
    plot.xaxis.axis_label_text_font_size = "12pt"
    plot.yaxis.axis_label = "UMAP 2"
    plot.yaxis.axis_label_text_font_size = "12pt"

    
    plot.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    plot.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
    plot.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    plot.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
    plot.xaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    plot.yaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels

    plot.add_layout(Title(text="Z-score", align="left",standoff=0,offset=110), "right")

    # Save static images 
    for static_image_i,static_image in enumerate(static_images_save):
        plot_scatter(x=umap_df['x'],y=umap_df['y'],values=values_dict[static_image],query=highlight_query,title=static_image,min_val=min_val,max_val=max_val,arrow_loc=umap_df.loc[highlight_query],rank=static_image_i+1,filename=file_path+'static/')
 

    default_text = "Figure {}. {}{}"
    pre = Paragraph(text = default_text.format(figure_counter, caption_text, first_selection), width=1000, height=10, style={"font-family":'Helvetica', "font-style": "italic"})
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
                source.change.emit();
            """)
            #plot.height = 1000+10*(category_list_dict[val].length)
        else:
            if len(option_list) >0:
                callback_adt = CustomJS(args=dict(source=source, \
                                                pre=pre, \
                                                values_dict=values_dict, \
                                                additional_info=additional_info,\
                                                figure_counter=figure_counter,
                                                color_mapper=color_mapper,\
                                                colors_dict=colors_dict,\
                                                scatter=scatter,
                                                caption_text=caption_text), code="""        
                    const val = cb_obj.value;    
                    source.data.values = values_dict[val]
                    if (additional_info != null) {
                        source.data.info = additional_info[val]
                    }
                    
                    pre.text = "Figure "+figure_counter+". "+caption_text+val+".";  
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
        select = Select(title="Color by " + color_by_title + ':', value=first_selection, options=option_list)
        select.js_on_change('value', callback_adt)
        
        col = column(select, row(column(plot, pre)))
        show(col)

    else:
        col = column(plot, pre)
        show(col)
    # Save to html file
    # output_file(filename=file_path+highlight_query+'_dynamic_umap.html', title=highlight_query)
    # save(col)

# Interactive network visualization
def network_vis(QUERY,LNCRNA_COEXP,GENES_2_ENSEMBL,ROW_GENES,NETWORK_EDGE_FILE,CHROM_LOC_DATA):
    
    edge_matrix = load_npz(NETWORK_EDGE_FILE)

    # Import chromosome location metadata
    chr_loc = pd.read_csv(CHROM_LOC_DATA,sep='\t')
    chr_loc['Chromosome/scaffold name'] = chr_loc['Chromosome/scaffold name'].astype(str)
    ensembl_2_chromsome = dict(zip(chr_loc['Gene stable ID'],chr_loc['Chromosome/scaffold name']))

    # Find chromosome locations for each node
    ensembl_network = [GENES_2_ENSEMBL[x] for x in LNCRNA_COEXP[0:100].index]
    ensembl_network_loc = []
    for x in ensembl_network:
        try:
            ensembl_network_loc.append(ensembl_2_chromsome[x])
        except:
            ensembl_network_loc.append('NA')
    
    # Node metadata
    network = LNCRNA_COEXP[0:100]
    #network.insert(1, 'Ensembl', ensembl_network, False)
    network.insert(1, 'Chromosome', ensembl_network_loc,False)

    # Load pre-computed edges
    # Correlations <0.3 are 0
   

    # Edges per node to start
    n_edges = 3

    # Find all pairwise correlations between top correlated genes
    network_genes = list(network.index)
    network_genes_done = copy.deepcopy(network_genes)
    idx_dict = dict(zip((ROW_GENES),list(range(0,len(ROW_GENES)))))

    node1 = []
    node2 = []
    edge_values = []
    
    for symbol1 in list(network_genes):
        idx1 = idx_dict[symbol1]
        for symbol2 in network_genes_done:
            idx2 = idx_dict[symbol2]
            c = edge_matrix[idx1,idx2]
            node1.append(symbol1)
            node2.append(symbol2)
            edge_values.append(c)
        network_genes_done.remove(symbol1)
        
    # Add the top n_edges for the lncRNA of interest - regardless of correlation
    for i in range(0,n_edges):
        node1.append(list(network.index)[i])
        node2.append(QUERY)
        edge_values.append(list(network["Pearson's Correlation Coefficient"])[i])

    # Edge Dataframe
    all_edges_df = pd.DataFrame({"Node1":node1,'Node2':node2,'Correlation':edge_values})
    all_edges_df = all_edges_df[all_edges_df['Correlation']!=1].reset_index(drop=True) # remove self-correlations
    all_edges_df = all_edges_df[all_edges_df['Correlation']!=0].reset_index(drop=True) # remove correlations of 0
    
    # Only keep the top n edges per node to start 
    idx_keep = []
    for symbol in list(network_genes)+list([QUERY]):
        top_idx=list(all_edges_df[(all_edges_df['Node1']==symbol) | (all_edges_df['Node2']==symbol)].sort_values(by='Correlation',ascending=False)[0:n_edges].index)
        for idx in top_idx:
            idx_keep.append(idx)
    all_edges_df = all_edges_df.loc[sorted(np.unique(idx_keep))]
    
    
    # Calculate the average edges per node in the network 
    average = []
    for node in np.unique(list(all_edges_df['Node1'])+list(all_edges_df['Node2'])):
        sub_df = all_edges_df[(all_edges_df['Node1']==node) | (all_edges_df['Node2']==node)]
        average.append(len(sub_df))
    average = np.mean(average)

    # Prune the network
    # Remove edges from hub node one by one until the network has an average <3 edges per node
    # Hub nodes are nodes with >10 edges. If there are no more hub nodes and the average is still >3 edges/node, decrease the definition of a hub node by 1
    hub_node_n = 10
    while average >3:
        n_hubs = 0
        remove_edges = []
        for symbol in list(network_genes):
            if symbol in np.unique(list(all_edges_df.Node1) + list(all_edges_df.Node2)):
                sub_hub = all_edges_df[(all_edges_df.Node1 == symbol) | (all_edges_df.Node2 == symbol)].sort_values(by='Correlation',ascending=False)
                if len(sub_hub)>hub_node_n:
                    n_hubs = n_hubs + 1
                    index_remove = list(sub_hub.index[hub_node_n-1:])
                    for idx in index_remove:
                        remove_edges.append(idx)
        remove_edges = np.unique(remove_edges)
        edges_keep = [x for x in list(all_edges_df.index) if x not in remove_edges]
        all_edges_df= all_edges_df.loc[edges_keep]
        average_new = []
        for node in np.unique(list(all_edges_df['Node1'])+list(all_edges_df['Node2'])):
            sub_df = all_edges_df[(all_edges_df['Node1']==node) | (all_edges_df['Node2']==node)]
            average_new.append(len(sub_df))
        average = np.mean(average_new)
        if n_hubs == 0:
            hub_node_n = hub_node_n-1
    all_edges_df= all_edges_df.reset_index(drop=True)
    
    # After pruning, add the top 5 edges for the lncRNA of interest regardless of the pruning step
    lncRNA_connections = all_edges_df[(all_edges_df['Node1']==QUERY)|(all_edges_df['Node2']==QUERY)]
    for i,node in enumerate(LNCRNA_COEXP[0:5].index):
        if node not in list(lncRNA_connections['Node1']):
            if node not in list(lncRNA_connections['Node2']):
                    df_add = {'Node1':node,'Node2':QUERY,'Correlation':LNCRNA_COEXP["Pearson's Correlation Coefficient"][i]}
                    all_edges_df = all_edges_df.append(df_add, ignore_index=True)


    # If a node has no edges after pruning, remove it
    nodes_keep = np.unique(list(all_edges_df['Node1'])+list(all_edges_df['Node2']))
    nodes_keep = [x for x in nodes_keep if x!=QUERY]
    network = network.loc[nodes_keep]
    
    # Get hover titles for nodes
    titles = []
    for i,x in enumerate(network.index):
        titles.append(str(x)+ "(Chromosome:" + network['Chromosome'][i]+")")
    titles.append(QUERY+ "(Chromosome:" + ensembl_2_chromsome[GENES_2_ENSEMBL[QUERY]]+")")

    # Get node colors
    color_list = list(Category20[20] + Accent[8][5:6]+ Category20b[20][0:1] + Category20b[20][4:5]+ Bokeh[8][7:]+Colorblind[8][5:6])
    color_list = [x+'90' for x in color_list] # Make node colors more opaque 
    color_palette = list(random.sample(list(itertools.chain(*zip(color_list))),len(list(np.unique(list(network['Chromosome'])+list(ensembl_2_chromsome[GENES_2_ENSEMBL[QUERY]]))))))
    chrom_2_color= dict(zip(list(np.unique(list(network['Chromosome'])+list(ensembl_2_chromsome[GENES_2_ENSEMBL[QUERY]]))),color_palette ))
    colors = [chrom_2_color[x] for x in list(network['Chromosome'])]
    colors.append('red') # The lncRNA is labeled red

    # Create the network
    g = Network(height=800,width=800, notebook=True,bgcolor='wh', font_color='black')
    g.add_nodes(
        list(list(network.index)+list([QUERY])),
        value=[10 for x in titles],
        title=titles,
        label=list(list(network.index)+list([QUERY])),
        color=colors
    )

    # Add edges to the network
    for i,node in enumerate(list(all_edges_df['Node1'])):
        g.add_edge(node, list(all_edges_df['Node2'])[i], value=list(all_edges_df['Correlation'])[i])

    # Set network options 
    g.set_options('''
        var options = {
        "nodes": {
            "color":{
            "border": "black"
            },
        "borderWidth": 3,
            "font": {
            "color":"black",
            "size": 20,
            "face": "arial",
            "highlight":"orange"
            },
            "shadow": {
            "enabled": true
            }
        },
        "edges": {
            "color": {
            "color":"lightgrey",
            "inherit": false,
            "highlight":"orange"
            },
            "smooth": false
        },
        "physics": {
            "hierarchicalRepulsion": {
            "centralGravity": 0
            },
            "minVelocity": 0.75
        }
        }
        ''')
    # Display network
    return(g,network,all_edges_df)
        