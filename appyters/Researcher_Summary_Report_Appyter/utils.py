#import packages
from collections import defaultdict
from IPython.display import display,FileLink, HTML, Markdown, IFrame
import requests
import plotly.express as px
import plotly
from scholarly import scholarly
import datetime
import math
from bs4 import BeautifulSoup
from PIL import Image, ImageDraw, ImageFont
import textwrap
import matplotlib.pyplot as plt
plotly.offline.init_notebook_mode()

# Class definition to define HTML object for the text display cards within the appyter.
class MyHTML:
  def __init__(self, html):
    self._html = html
  def _repr_html_(self):
    return self._html
class MyMarkdown:
  def __init__(self, markdown):
    self._markdown = markdown
  def _repr_markdown_(self):
    return self._markdown
#255,255,255
def make_bar_plot(input_data_dict, x_axis_title, y_axis_title, graph_title, source_title=''):
    graph_title = '<br>'.join(textwrap.wrap(graph_title, 40))
    dict_for_graph = defaultdict(list)
    for key, value in input_data_dict.items():
        dict_for_graph[x_axis_title].append(key)
        dict_for_graph[y_axis_title].append(value)
    fig = px.bar(dict_for_graph, x = x_axis_title, y=y_axis_title, title=graph_title)
    update_layout_params = dict(
        font = dict(size = 12),
        plot_bgcolor="rgba(240, 238, 240,0)",
        title_font_size=32,
        yaxis=dict(title_font=dict(size=18)), 
        width = 750, 
        height = 500)
    #Making sure the the x axis is linear when we have less than 35 data points
    if len(input_data_dict.keys()) < 35:
        fig.update_layout(xaxis = dict(
        tickmode = 'linear',
        dtick = 1, 
        title_font=dict(size=18)
        ))
    #Otherwise skipping every 2 so the x axis is more readable. 
    else:
        fig.update_layout(xaxis = dict(
        tickmode = 'linear',
        dtick = 2, 
        title_font=dict(size=18)
        ))
    fig.update_layout(update_layout_params)
    fig.update_traces(marker_color='black')
    fig.update_traces(width=0.99)
    #Adding the source for all the plots except the drugshot and geneshot ones. 
    if type(list(input_data_dict.keys())[0]) != str:
        fig.update_layout(xaxis=dict(range=[min(input_data_dict.keys())-0.5, max(input_data_dict.keys())+0.5]))
            # add annotation
        #Annotation at bottom right
        fig.add_annotation(dict(font=dict(color='black',size=12),
                                        x= 1.1,
                                        y = -.15,
                                        showarrow=False,
                                        text=source_title,
                                        textangle=0,
                                        xanchor='right',
                                        yanchor='top',
                                        xref="paper",
                                        yref="paper")) 
    else:
        #Annotation on top right
        fig.add_annotation(dict(font=dict(color='black',size=12),
                                        x= 1.1,
                                        y = 1.32,
                                        showarrow=False,
                                        text=source_title,
                                        textangle=0,
                                        xanchor='right',
                                        yanchor='top',
                                        xref="paper",
                                        yref="paper")) 
    #Fixing the y axis of the graphs when the values are on a smaller scale but still greater than frequency between 0 and 1. 
    if max(dict_for_graph[y_axis_title]) <= 10 and max(dict_for_graph[y_axis_title]) > 1:
        fig.update_layout(yaxis=dict(dtick=1))
    return fig

def make_line_plot(input_data_dict, x_axis_title, y_axis_title, graph_title, source_title=''):
    graph_title = '<br>'.join(textwrap.wrap(graph_title, 40))
    dict_for_graph = defaultdict(list)
    running_sum = 0
    for key, value in input_data_dict.items():
        running_sum += value
        dict_for_graph[x_axis_title].append(key)
        dict_for_graph[y_axis_title].append(running_sum)
    fig = px.line(dict_for_graph, x = x_axis_title, y=y_axis_title, title=graph_title, markers = True)
    update_layout_params = dict(
        font = dict(size = 12),
        plot_bgcolor="rgba(240, 238, 240,0)",
        title_font_size=32,
        yaxis=dict(title_font=dict(size=18)), 
        width = 750, 
        height = 500)
    fig.update_layout(update_layout_params)
    if len(input_data_dict.keys()) < 35:
        fig.update_layout(xaxis = dict(
        tickmode = 'linear',
        dtick = 1, 
        title_font=dict(size=18)
        ))
    else:
        fig.update_layout(xaxis = dict(
        tickmode = 'linear',
        dtick = 2, 
        title_font=dict(size=18)
        ))
    #Line color change
    fig.update_traces(line_color='black')
    if type(list(input_data_dict.keys())[0]) != str:
        fig.update_layout(xaxis=dict(range=[min(input_data_dict.keys())-0.5, max(input_data_dict.keys())+0.5]))
        fig.add_annotation(dict(font=dict(color='black',size=12),
                                x= 1.1,
                                y = -.15,
                                showarrow=False,
                                text=source_title,
                                textangle=0,
                                xanchor='right',
                                yanchor='top',
                                xref="paper",
                                yref="paper")) 
    if max(dict_for_graph[y_axis_title]) <= 10 and max(dict_for_graph[y_axis_title]) > 1:
        fig.update_layout(yaxis=dict(dtick=1))
    return fig
def calculate_ar_index(citations):
    """
    Calculate the AR index as based in the paper from BiHui et al. 
    AR-index calculated by taking into account the years that passed
    """
    h_index = 0
    running_sum = 0
    year_now = int(datetime.date.today().year)
    for i, count in enumerate(citations):
        if count[2] >= i+1:
            h_index = i+1
            running_sum += count[2]/(year_now-count[1] + 1)
        else:
            break
    # print(f"The calculated h-index is: {h_index}")
    return math.sqrt(running_sum)

#Displaying the figure number and if a title is given, then displaying the download link. 
def display_figure_labels(output_folder, counter, caption, title = None):
    display(MyMarkdown("*Figure {}. {}*".format(counter, caption)))
    if title != None:
        display(FileLink(output_folder+title+'.png' , result_html_prefix=str('Download Figure {} (PNG): '.format(counter))))
    counter += 1
    return counter


def query_semantic_scholar_citation(name_of_researcher):
    citation_dict = defaultdict(int)
    #Retrieving the most likely author id based off the name of the input passed in from the semantic Scholar API. 
    url = f"https://api.semanticscholar.org/graph/v1/author/search?query={name_of_researcher}&fields=paperCount,name"
    researcher_query_response = requests.get(url)
    if researcher_query_response.status_code == 200:
        data = researcher_query_response.json()
        id_of_researcher = None
        final_name = None
        running_count = 0
        affil = None
        #Loop through authors returned from the API Call. 
        if data['total'] > 0:
            for dict_author in data['data']:
                id_now = dict_author['authorId']
                name_returned = dict_author['name']
                paper_count = dict_author['paperCount']
                name_returned = name_returned.replace("’", "'")
                #Currently using the author with the greatest number of papers and saving their information
                if paper_count > running_count and name_of_researcher.split()[-1] in name_returned:
                    id_of_researcher = id_now
                    running_count = paper_count
                    final_name = name_returned.replace(" ", "-")
            #Setting up API request to the author endpoint in order to get the citation counts for each paper and year of the researcher. 
            url_for_papers_final = f"https://api.semanticscholar.org/graph/v1/author/{id_of_researcher}?fields=name,citationCount,paperCount,hIndex,aliases,papers.year,papers.citationCount"
            citation_response = requests.get(url_for_papers_final)
            if citation_response.status_code == 200:
                display(MyMarkdown("### Link to [Semantic Scholar Page](https://www.semanticscholar.org/author/{}/{}) for {}".format(final_name, id_of_researcher, name_of_researcher)))
                data = citation_response.json()
                #If there are more than 500 papers, have to go to the papers endpoint with a lst of papers in order to get that information since we have to use offsets
                if data['paperCount'] > 500:
                    dict_for_holding_ids = {}
                    total = 0
                    offset = 0
                    key_for_dict = 0
                    while total < data['paperCount']:
                        url_for_papers_endpoint = "https://api.semanticscholar.org/graph/v1/author/{}/papers?fields=citations.year&offset={}&limit=500".format(id_of_researcher, str(offset))
                        res = requests.get(url_for_papers_endpoint)
                        paper_list_data = res.json()
                        string_of_ids = ""
                        array_holding_ids = []
                        if 'data' in paper_list_data:
                            for paper in paper_list_data['data']:
                                string_of_ids += paper['paperId'] + ","
                                array_holding_ids.append(paper['paperId'])
                            dict_for_holding_ids[key_for_dict] = array_holding_ids
                            key_for_dict += 1
                            total += len(paper_list_data['data'])
                            if 'next' in paper_list_data:
                                offset = paper_list_data['next']
                    for key_ in dict_for_holding_ids:
                        url_for_multiple_paper_info = 'https://api.semanticscholar.org/graph/v1/paper/batch'
                        parameters = {'fields': 'year,citationCount,title'}
                        data_dict = {"ids": dict_for_holding_ids[key_]}
                        response = requests.post(url_for_multiple_paper_info,params=parameters,json=data_dict)
                        output_data = response.json()
                        for item in output_data:
                            if 'year' and 'citationCount' in item:
                                if item['year'] != None:
                                    citation_dict[item['year']] += item['citationCount']
                #Add the year and citation count information if there are less than 500 papers. 
                else:
                    for paper in data['papers']:
                        if type(paper['year']) == int:
                            citation_dict[paper['year']] += paper['citationCount']
            #Error with the id of the reseacher in the API call. 
            else:
                citation_dict = None
                print("Error in getting the citation information for researcher by ID.")

        if citation_dict != None and len(citation_dict) != 0:
            # print(citation_dict)
            year_keys = list(citation_dict.keys())
            year_keys.sort()
            citation_dict = {year:citation_dict[year] for year in year_keys}
            return citation_dict
        return citation_dict
    else:
        print('Error in querying this researcher from Semantic Scholar by name. Their information may not be here. A manual search may help. ')
        return citation_dict


def getting_information_from_openalex(name_of_researcher,output_folder):
    url_link = "https://api.openalex.org/authors?search={}".format(name_of_researcher)
    response = requests.get(url_link)
    if response.status_code == 200:
        interests = []
        h_index = None
        i10_index = None
        total_times_cited = None
        institution = ''

        data = response.json()
        try:
            if data['meta']['count'] > 0:
                first_result = data['results'][0]
                if name_of_researcher.split()[0] in first_result['display_name'] and name_of_researcher.split()[-1] in first_result['display_name']:
                    total_times_cited = first_result['cited_by_count']
                    if 'h_index' in first_result['summary_stats']:
                        h_index = first_result['summary_stats']['h_index']
                    if 'i10_index' in first_result['summary_stats']:
                        i10_index = first_result['summary_stats']['i10_index']
                    if 'last_known_institution' in first_result and first_result['last_known_institution'] != None:
                        institution = first_result['last_known_institution']['display_name']
                    for concept in first_result['x_concepts'][:5]:
                        interests.append(concept['display_name'])
                else:
                    print("Name returned from OpenAlex API doesnt match name of researcher.")
                
                # display_summary_text_from_openalex(institution,interests, h_index,i10_index, total_times_cited, name_of_researcher)
                return display_summary_text_from_openalex_png(institution,interests, h_index,i10_index, total_times_cited, name_of_researcher,output_folder)

            else:
                print("No matching information for this researcher from OpenAlex")
        except:
            print("Error in the OpenAlex response JSON")
            return None
    else:
        print("Error in querying from OpenAlex API.")
    return None


def getting_information_from_wiki(name_of_researcher, output_folder):
    wiki_institution = ''
    wiki_known_for = []
    wiki_field_interests = []
    wiki_awards = []
    wiki_search_url = f"http://en.wikipedia.org/w/api.php?action=query&list=search&srsearch={name_of_researcher}&format=json"
    title_to_get_content_from = None
    page_id = None
    wiki_title_response = requests.get(wiki_search_url)
    if wiki_title_response.status_code == 200:
        data = wiki_title_response.json()
        if data['query']['searchinfo']['totalhits'] > 0:
            for page_query in data['query']['search']:
                if name_of_researcher.split()[0] in page_query['title'] and name_of_researcher.split()[-1] in page_query['title']:
                    title_to_get_content_from = page_query['title']
                    page_id = page_query['pageid']
                    break
            if title_to_get_content_from != None:
                display(MyMarkdown("### Link to [Wiki Page](https://en.wikipedia.org/wiki/{}) for {}".format(title_to_get_content_from.replace(" ","_"), name_of_researcher)))
                url_to_get_page_content = f"https://en.wikipedia.org/w/api.php?action=parse&page={title_to_get_content_from}&props=text&format=json"
                response = requests.get(url_to_get_page_content)
                data  = response.json()
                try:
                    soup = BeautifulSoup(data['parse']['text']['*'], 'html.parser')
                    table_data = soup.find('table')
                    if table_data != None:
                        table_row_data = table_data.find_all('tr')
                        for data in table_row_data:
                            th = data.find('th')
                            td = data.find('td')
                            if td != None and th != None:
                                if "Known" in th.text:
                                    for content in td.contents:
                                        if type(content) == type(BeautifulSoup("", "html.parser").new_tag("tag")) and len(content.text) > 0:
                                            # print(content.text)
                                            wiki_known_for.append(content.text)
                                elif "Fields" in th.text:
                                    for content in td.contents:
                                        if type(content) == type(BeautifulSoup("", "html.parser").new_tag("tag")) and len(content.text) > 0:
                                            wiki_field_interests.append(content.text)
                                            # print(content.text)
                                elif "Institutions" in th.text:
                                    for index, content in enumerate(td.contents):
                                        if type(content) == type(BeautifulSoup("", "html.parser").new_tag("tag")) and len(content.text) > 0:
                                            # print(content.text)
                                            if index != len(td.contents) - 1:
                                                wiki_institution += content.text + " and "
                                            else:
                                                wiki_institution += content.text
                                elif "Awards" in th.text:
                                    for content in td.contents:
                                        if type(content) == type(BeautifulSoup("", "html.parser").new_tag("tag")) and len(content.text) > 0:
                                            wiki_awards.append(content.text)
                                            # print(content.text)
                        # display_summary_text_from_wikipedia(wiki_institution, wiki_known_for, wiki_field_interests, wiki_awards, name_of_researcher)
                        return display_summary_text_from_wikipedia_png(wiki_institution, wiki_known_for, wiki_field_interests, wiki_awards, name_of_researcher, output_folder)
                    else:
                        print("No data to show")
                except:
                    print("Error in reading the content from Wikipedia page")
                    return None
            else:
                print("No matching Wikipedia Page for this researcher")
        else:
            print("No matching Wikipedia Page for this researcher")
    else:
        print("Error in querying Wikipedia for the name of researcher")   
    return None      
#Old image color (221, 224, 237)
#Making the png image cards with text displayed on top of the png cards and saving those pngs as output. 
font_for_title = 'LiberationSans-Bold.ttf'
font_for_content = 'LiberationSans-Regular.ttf'
# font_for_title = 'OpenSans-Bold.ttf'
# font_for_content = 'OpenSans-Regular.ttf'
title_font_size_ = 18
content_font_size_ = 16
source_font_size_ = 10
image_width_  = 200
image_height_ = 200
dpi_value_ = (200,200)
def display_summary_text_from_wikipedia_png(wiki_institution = '', wiki_known_for = [], wiki_field_interests = [], wiki_awards = [], name_of_researcher = None, output_folder = None):
    display(MyMarkdown("## Text Summary Information for {} (From Wikipedia) ##".format(name_of_researcher)))
    display_list = []
    title_font_size = title_font_size_
    content_font_size = content_font_size_
    image_width = image_width_
    image_height = image_height_
    title_font = ImageFont.truetype(font_for_title, title_font_size)
    content_font = ImageFont.truetype(font_for_content, content_font_size)
    source_font_size = source_font_size_
    source_font = ImageFont.truetype(font_for_content, source_font_size)
    source_title = 'Sourced from Wikipedia'
    dpi_value = dpi_value_
    if wiki_institution != '':
        inst_image = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(inst_image)
        text_1 = textwrap.wrap("Institution Affiliation", width = 25)
        text_2 = textwrap.wrap(wiki_institution, width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size

        for line in text_2:
            draw.text((5, y+10), line, font=content_font, fill='black')
            y += content_font_size
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        inst_image.save(output_folder + "1card_wiki_image_1.png", dpi=dpi_value)
        display_list.append(inst_image)
        # display(inst_image)
    if len(wiki_awards) > 0: 
        awards = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(awards)
        text_1 = textwrap.wrap("Awards", width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font= title_font, fill='black')
            y += content_font_size
        for award in wiki_awards:
            text_2 = textwrap.wrap(award, width = 20)
            for line in text_2:
                draw.text((5, y+10), line, font = content_font, fill='black')
                y += content_font_size
            y += 10
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        awards.save(output_folder + "1card_wiki_image_2.png", dpi=dpi_value)
        display_list.append(awards)
        # display(awards)

    if len(wiki_known_for) > 0:
        known_for = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(known_for)
        text_1 = textwrap.wrap("Known For", width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font=title_font, fill='black')
            y += content_font_size
        for known in wiki_known_for:
            text_2 = textwrap.wrap(known, width = 20)
            for line in text_2:
                draw.text((5, y+10), line, font=content_font, fill='black')
                y += content_font_size
            y += 10
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        known_for.save(output_folder + "1card_wiki_image_3.png", dpi=dpi_value)
        display_list.append(known_for)
        # display(known_for)

    if len(wiki_field_interests) > 0:
        interests = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(interests)
        text_1 = textwrap.wrap("Research Fields", width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font=title_font, fill='black')
            y += content_font_size
        for interest in wiki_field_interests:
            text_2 = textwrap.wrap(interest, width = 20)
            for line in text_2:
                draw.text((5, y+10), line, font=content_font, fill='black')
                y += content_font_size
            y += 10
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        interests.save(output_folder + "1card_wiki_image_4.png", dpi=dpi_value)
        display_list.append(interests)
        # display(interests)
    #Creating a 1x3 figure with 1-based indexing to add to the add_subplot function
    fig = plt.figure(figsize=(10, 10))
    rows = 1
    columns = len(display_list)
    for idx, img in enumerate(display_list):
        fig.add_subplot(rows, columns, idx + 1)
        plt.imshow(img)
        plt.axis('off')
    
    plt.show()
    return display_list


def display_summary_text_from_openalex_png(institution = '', interests = [], h_index = None, i10_index = None, total_times_cited = None, name_of_researcher = None, output_folder = None):
    display_list = []
    title_font_size = title_font_size_
    content_font_size = content_font_size_
    image_width = image_width_
    image_height = image_height_
    title_font = ImageFont.truetype(font_for_title, title_font_size)
    content_font = ImageFont.truetype(font_for_content, content_font_size)
    source_font_size = source_font_size_
    source_font = ImageFont.truetype(font_for_content, source_font_size)
    source_title = 'Sourced from OpenAlex'
    dpi_value = dpi_value_
    display(MyMarkdown("## Text Summary Information for {} (From OpenAlex API) ##".format(name_of_researcher)))
    if institution != '':
        inst_image = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(inst_image)
        text_1 = textwrap.wrap("Institution Affiliation", width = 25)
        text_2 = textwrap.wrap(institution, width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size

        for line in text_2:
            draw.text((5, y+10), line, font=content_font, fill='black')
            y += content_font_size
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        inst_image.save(output_folder + "1card_openalex_image_1.png", dpi=dpi_value)
        display_list.append(inst_image)
        # display(inst_image)
    if h_index != None:
        h_image = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(h_image)
        title = 'Metrics'
        h_index = textwrap.wrap(f'H-index: {h_index}', width = 25)
        i10_index = textwrap.wrap(f'I10-index : {i10_index}', width = 25)
        times_cited = textwrap.wrap('Total Citations: '+str(total_times_cited), width=25)
        y = 0
        draw.text((5, y), title, font=title_font, fill='black')
        y += title_font_size + 10
        for line in h_index:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        y += 10
        for line in i10_index:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        y += 10
        for line in times_cited:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        h_image.save(output_folder + "1card_openalex_image_2.png", dpi=dpi_value)
        display_list.append(h_image)
        # display(h_image)
    if len(interests) > 0:
        interests_image = Image.new("RGB", (image_width, image_height), (240, 238, 240))
        draw = ImageDraw.Draw(interests_image)
        title = textwrap.wrap('Research Interests', width = 25)
        y = 0
        for line in title:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size
        y += 15
        for interest in interests:
            interest_text = textwrap.wrap(interest, width=20)
            for line in interest_text:
                draw.text((5, y), line, font=content_font, fill='black')
                y += content_font_size
            y += 10
        draw.text((image_width*0.35, image_height-15), source_title, font=source_font, fill='black')
        interests_image.save(output_folder + "1card_openalex_image_3.png", dpi=dpi_value)
        display_list.append(interests_image)
        # display(interests_image)
    
    fig = plt.figure(figsize=(10, 10))
    rows = 1
    columns = len(display_list)
    for idx, img in enumerate(display_list):
        fig.add_subplot(rows, columns, idx + 1)
        plt.imshow(img)
        plt.axis('off')
    plt.show()
    return display_list

