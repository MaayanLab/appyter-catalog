
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

def plot_bar(df,title,x_label,y_label):
    df_dict = {'x':list(df.index),'y':df[list(df.columns)[0]]}
    bar = figure(x_range=df_dict['x'], height=500, width=1000, title=title,toolbar_location="right", tools=["hover",'save'],tooltips=[ (x_label, "@x"),(y_label, "@y")])
    bar.vbar(x='x', top='y', width=0.7, source=df_dict)
    bar.xgrid.grid_line_color = None
    #bar.y_range.start = 0
    bar.xaxis.major_label_orientation = 1
    bar.xaxis.axis_label = x_label
    bar.yaxis.axis_label = y_label
    show(bar)

# Plot the top terms for each prediction library
def plot_results(library_names, results_dfs, top_results=20):
    
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
        text = go.Scatter(
            x=[max(bar['x'])/50 for x in range(len(bar['y']))],
            y=bar['y'],
            mode='text',
            hoverinfo='none',
            showlegend=False,
            text=['<b>{}</b>'.format(rowData['Term']) for index, rowData in results_df[0:top_results].iterrows()],
            textposition="middle right",
            textfont={'color': 'black','size':8})
        fig.append_trace(text, 1, i+1)
    
    annotations= [{'x': 0.25, 'y': 1.1, 'text': '<span style="color: black; font-size: 15pt; font-weight: 600;">' +library_names[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},{'x': 0.75, 'y': 1.1, 'text': '<span style="color: black; font-size: 15pt; font-weight: 600;">' +library_names[1]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}]
    fig['layout'].update(height = 500, hovermode='closest', annotations=annotations)
    fig.update_layout(title='',height = 500,title_font_size = 25,title_x=0.5)
    
    fig['layout']['xaxis1'].update(domain=[0, 0.49], title='Mean Pearson Correlation' ,range=(0,max_scores[0]+max_scores[0]*.01))
    fig['layout']['xaxis2'].update(domain=[0.51, 1], title='Mean Pearson Correlation',range=(0,max_scores[1]+max_scores[1]*.01))
    fig['layout']['yaxis1'].update(showticklabels=False)
    fig['layout']['yaxis2'].update(showticklabels=False)
    fig['layout']['margin'].update(l=30, t=65, r=30, b=35)
    
    iplot(fig)

def str_to_int(string, mod):
    string = re.sub(r"\([^()]*\)", "", string).strip()
    byte_string = bytearray(string, "utf8")
    return int(hashlib.sha256(byte_string).hexdigest(), base=16)%mod

def plot_scatter(umap_df, values_dict, option_list, sample_names, caption_text, category_list_dict=None, location='right', category=True, dropdown=False, figure_counter=0, additional_info=None, color_by_title='', first_selection = None, highlight_query=None):
    if first_selection == None:
        first_selection = option_list[0]
    # init plot 
    source = ColumnDataSource(data=dict(x=umap_df["x"], y=umap_df["y"], values=values_dict[first_selection], names=sample_names))
    
    # node size
    if umap_df.shape[0] > 1000:
        node_size = 4
    else:
        node_size = 6
    plot = figure(plot_width=1000, plot_height=800)  
 
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

        color_mapper = LinearColorMapper(palette=colors_dict[first_selection] , low=min(values_dict[first_selection]), high=max(values_dict[first_selection]))
        color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12)
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
    default_text = "Figure {}. {}{}"
    pre = Paragraph(text = default_text.format(figure_counter, caption_text, first_selection), width=1000, height=20, style={"font-family":'Helvetica', "font-style": "italic"})
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