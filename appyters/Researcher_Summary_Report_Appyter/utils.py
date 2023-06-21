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
    
def make_bar_plot(input_data_dict, x_axis_title, y_axis_title, graph_title, source_title=''):
    dict_for_graph = defaultdict(list)
    for key, value in input_data_dict.items():
        dict_for_graph[x_axis_title].append(key)
        dict_for_graph[y_axis_title].append(value)
    fig = px.bar(dict_for_graph, x = x_axis_title, y=y_axis_title, title=graph_title)
    update_layout_params = dict(
        font = dict(size = 12),
        plot_bgcolor="rgba(255,255,255,0)",
        title_font_size=18,
        yaxis=dict(title_font=dict(size=18)), 
        width = 750, 
        height = 500)
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
    fig.update_layout(update_layout_params)
    fig.update_traces(marker_color='black')
    fig.update_traces(width=0.99) 
    if type(list(input_data_dict.keys())[0]) != str:
        fig.update_layout(xaxis=dict(range=[min(input_data_dict.keys())-0.5, max(input_data_dict.keys())+0.5]))
            # add annotation
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

def make_line_plot(input_data_dict, x_axis_title, y_axis_title, graph_title, source_title=''):
    dict_for_graph = defaultdict(list)
    running_sum = 0
    for key, value in input_data_dict.items():
        running_sum += value
        dict_for_graph[x_axis_title].append(key)
        dict_for_graph[y_axis_title].append(running_sum)
    fig = px.line(dict_for_graph, x = x_axis_title, y=y_axis_title, title=graph_title, markers = True)
    update_layout_params = dict(
        font = dict(size = 12),
        plot_bgcolor="rgba(255,255,255,0)",
        title_font_size=18,
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
def display_figure_labels(counter, caption, title = None):
    display(MyMarkdown("*Figure {}. {}*".format(counter, caption)))
    if title != None:
        display(FileLink("output_images/"+title+'.png' , result_html_prefix=str('Download Figure {} (PNG): '.format(counter))))
    counter += 1
    return counter


def query_google_citation(name_of_researcher):
    search_query = scholarly.search_author(name_of_researcher)
    try:
        author  = next(search_query)
        display(MyMarkdown("### Link to [Google Scholar Page](https://scholar.google.com/citations?user={}) for {}".format(author['scholar_id'], name_of_researcher)))
        summary_info = scholarly.fill(author, sections=['counts', 'indices'])
        # print(summary_info)
        author = scholarly.fill(author)
        # print(len([pub['bib']['title'] for pub in author['publications']]))
        list_storing_citations_and_years = []
        for pub in author['publications']:
            # print(pub)
            if 'pub_year'in pub['bib']:
                list_storing_citations_and_years.append([pub['bib']['title'], int(pub['bib']['pub_year']), int(pub['num_citations'])])
        # print(len(list_storing_citations_and_years))
        citation_dict = summary_info['cites_per_year']
        affiliation_from_google_scholar = summary_info['affiliation']
        interests_from_google_scholar = summary_info['interests']
        total_times_cited = summary_info['citedby']
        h_index_from_google_scholar = summary_info['hindex']
        h_index_from_google_scholar_last_5 = summary_info['hindex5y']

        list_storing_citations_and_years = sorted(list_storing_citations_and_years, key = lambda x:(x[2], x[1]), reverse=True)
        ar_index = calculate_ar_index(list_storing_citations_and_years)
        # display_summary_text_from_google_scholar(affiliation_from_google_scholar, h_index_from_google_scholar, interests_from_google_scholar, h_index_from_google_scholar_last_5, ar_index, total_times_cited, name_of_researcher)
        display_list = display_summary_text_from_google_scholar_png(affiliation_from_google_scholar, h_index_from_google_scholar, interests_from_google_scholar, h_index_from_google_scholar_last_5, ar_index, total_times_cited, name_of_researcher)
        return (citation_dict, display_list)
    except:
        print("No Google Scholar information for {}".format(name_of_researcher))
        return None


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
        if data['total'] > 0:
            for dict_author in data['data']:
                id_now = dict_author['authorId']
                name_returned = dict_author['name']
                paper_count = dict_author['paperCount']
                name_returned = name_returned.replace("’", "'")
                if paper_count > running_count and name_of_researcher.split()[-1] in name_returned:
                    id_of_researcher = id_now
                    running_count = paper_count
                    final_name = name_returned.replace(" ", "-")
            url_for_papers_final = f"https://api.semanticscholar.org/graph/v1/author/{id_of_researcher}?fields=name,citationCount,paperCount,hIndex,aliases,papers.year,papers.citationCount"
            citation_response = requests.get(url_for_papers_final)
            if citation_response.status_code == 200:
                display(MyMarkdown("### Link to [Semantic Scholar Page](https://www.semanticscholar.org/author/{}/{}) for {}".format(final_name, id_of_researcher, name_of_researcher)))
                data = citation_response.json()
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

                else:
                    for paper in data['papers']:
                        if type(paper['year']) == int:
                            citation_dict[paper['year']] += paper['citationCount']
            else:
                citation_dict = None
                print("Error in getting the citation information for researcher.")

        if citation_dict != None and len(citation_dict) != 0:
            # print(citation_dict)
            year_keys = list(citation_dict.keys())
            year_keys.sort()
            citation_dict = {year:citation_dict[year] for year in year_keys}
            #Make the citation bar and line graph
            return citation_dict
        return citation_dict
    else:
        print('Error in querying this researcher from Semantic Scholar. Their information may not be here. A manual search may help. ')
        return citation_dict


def getting_information_from_openalex(name_of_researcher):
    url_link = "https://api.openalex.org/authors?search={}".format(name_of_researcher)
    response  = requests.get(url_link)
    if response.status_code == 200:
        interests = []
        h_index = None
        i10_index = None
        total_times_cited = None
        institution = ''

        data = response.json()
        if data['meta']['count'] > 0:
            first_result = data['results'][0]
            if name_of_researcher.split()[0] in first_result['display_name'] and name_of_researcher.split()[-1] in first_result['display_name']:
                total_times_cited = first_result['cited_by_count']
                if 'h_index' in first_result['summary_stats']:
                    h_index = first_result['summary_stats']['h_index']
                if 'i10_index' in first_result['summary_stats']:
                    i10_index = first_result['summary_stats']['i10_index']
                if 'last_known_institution' in first_result:
                    institution = first_result['last_known_institution']['display_name']
                for concept in first_result['x_concepts'][:5]:
                    interests.append(concept['display_name'])
            
            # display_summary_text_from_openalex(institution,interests, h_index,i10_index, total_times_cited, name_of_researcher)
            return display_summary_text_from_openalex_png(institution,interests, h_index,i10_index, total_times_cited, name_of_researcher)

        else:
            print("No matching information for this researcher from OpenAlex")
    else:
        print("Error in querying from OpenAlex API.")
    return None


def getting_information_from_wiki(name_of_researcher):
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
                        return display_summary_text_from_wikipedia_png(wiki_institution, wiki_known_for, wiki_field_interests, wiki_awards, name_of_researcher)
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


def display_summary_text_from_wikipedia_png(wiki_institution = '', wiki_known_for = [], wiki_field_interests = [], wiki_awards = [], name_of_researcher = None):
    display(MyMarkdown("## Text Summary Information for {} (From Wikipedia) ##".format(name_of_researcher)))
    display_list = []
    title_font_size = 18
    content_font_size = 16
    image_width = 200
    image_height = 200
    title_font = ImageFont.truetype("Arial.ttf", title_font_size)
    content_font = ImageFont.truetype("Arial.ttf", content_font_size)
    if wiki_institution != '':
        inst_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
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
        inst_image.save("./output_images/wiki_image_1.png")
        display_list.append(inst_image)
        # display(inst_image)
    if len(wiki_awards) > 0: 
        awards = Image.new("RGB", (image_width, image_height), (221, 224, 237))
        draw = ImageDraw.Draw(awards)
        text_1 = textwrap.wrap("Awards", width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font= title_font, fill='black')
            y += content_font_size
        for award in wiki_awards:
            text_2 = textwrap.wrap(award, width = 20)
            for line in text_2:
                draw.text((5, y+10), line, font= content_font, fill='black')
                y += content_font_size
            y += 10
        awards.save("./output_images/wiki_image_2.png")
        display_list.append(awards)
        # display(awards)

    if len(wiki_known_for) > 0:
        known_for = Image.new("RGB", (image_width, image_height), (221, 224, 237))
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
        known_for.save("./output_images/wiki_image_3.png")
        display_list.append(known_for)
        # display(known_for)

    if len(wiki_field_interests) > 0:
        interests = Image.new("RGB", (image_width, image_height), (221, 224, 237))
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
        interests.save("./output_images/wiki_image_4.png")
        display_list.append(interests)
        # display(interests)

    fig = plt.figure(figsize=(10, 10))
    rows = 1
    columns = len(display_list)
    for idx, img in enumerate(display_list):
        fig.add_subplot(rows, columns, idx + 1)
        plt.imshow(img)
        plt.axis('off')
    
    plt.show()
    return display_list



def display_summary_text_from_google_scholar_png(affiliation_from_google_scholar = None, h_index_from_google_scholar = None, interests_from_google_scholar = [], h_index_from_google_scholar_last_5 = None, ar_index = None, total_times_cited = None, name_of_researcher = None):
    display(MyMarkdown("## Text Summary Information for {} (From Google Scholar) ##".format(name_of_researcher)))
    display_list = []
    title_font_size = 18
    content_font_size = 16
    image_width = 200
    image_height = 200
    title_font = ImageFont.truetype("Arial.ttf", title_font_size)
    content_font = ImageFont.truetype("Arial.ttf", content_font_size)
    if affiliation_from_google_scholar != '':
        inst_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
        draw = ImageDraw.Draw(inst_image)
        text_1 = textwrap.wrap("Institution Affiliation", width = 25)
        text_2 = textwrap.wrap(affiliation_from_google_scholar, width = 20)
        y = 0
        for line in text_1:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size

        for line in text_2:
            draw.text((5, y+10), line, font=content_font, fill='black')
            y += content_font_size
        inst_image.save("./output_images/googlescholar_image_1.png")
        display_list.append(inst_image)
        # display(inst_image)
    if h_index_from_google_scholar != None:
        h_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
        draw = ImageDraw.Draw(h_image)
        title = 'Metrics'
        h_index = textwrap.wrap(f'H-index: {h_index_from_google_scholar}', width = 25)
        h_index_5 = textwrap.wrap(f'H-index over last 5 years: {h_index_from_google_scholar_last_5}', width = 25)
        ar_index = textwrap.wrap(f'Age Related Index: {int(ar_index)}', width=25)
        times_cited = textwrap.wrap('Total Citations: '+str(total_times_cited), width=25)
        y = 0
        draw.text((5, y), title, font=title_font, fill='black')
        y += title_font_size + 10
        for line in h_index:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        y += 10
        for line in h_index_5:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        y += 10
        for line in ar_index:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        y += 10
        for line in times_cited:
            draw.text((5, y), line, font=content_font, fill='black')
            y += content_font_size
        h_image.save("./output_images/googlescholar_image_2.png")
        display_list.append(h_image)
        # display(h_image)
    if len(interests_from_google_scholar) > 0:
        interests_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
        draw = ImageDraw.Draw(interests_image)
        title = textwrap.wrap('Research Interests', width = 25)
        y = 0
        for line in title:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size
        y += 15
        for interest in interests_from_google_scholar:
            interest_text = textwrap.wrap(interest, width=25)
            for line in interest_text:
                draw.text((5, y), line, font=content_font, fill='black')
                y += content_font_size
            y += 10
        interests_image.save("./output_images/googlescholar_image_3.png")
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



def display_summary_text_from_openalex_png(institution = '', interests = [], h_index = None, i10_index = None, total_times_cited = None, name_of_researcher = None):
    display_list = []
    title_font_size = 18
    content_font_size = 16
    image_width = 200
    image_height = 200
    title_font = ImageFont.truetype("Arial.ttf", title_font_size)
    content_font = ImageFont.truetype("Arial.ttf", content_font_size)
    display(MyMarkdown("## Text Summary Information for {} (From OpenAlex API) ##".format(name_of_researcher)))
    if institution != '':
        inst_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
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
        inst_image.save("./output_images/openalex_image_1.png")
        display_list.append(inst_image)
        # display(inst_image)
    if h_index != None:
        h_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
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
        h_image.save("./output_images/openalex_image_2.png")
        display_list.append(h_image)
        # display(h_image)
    if len(interests) > 0:
        interests_image = Image.new("RGB", (image_width, image_height), (221, 224, 237))
        draw = ImageDraw.Draw(interests_image)
        title = textwrap.wrap('Research Interests', width = 25)
        y = 0
        for line in title:
            draw.text((5, y), line, font=title_font, fill='black')
            y += title_font_size
        y += 15
        for interest in interests:
            interest_text = textwrap.wrap(interest, width=25)
            for line in interest_text:
                draw.text((5, y), line, font=content_font, fill='black')
                y += content_font_size
            y += 10
        interests_image.save("./output_images/openalex_image_3.png")
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

def display_matrix(google_scholar_info, wiki_info, open_alex_info):
    wiki_info = wiki_info[:3]
    fig = plt.figure(figsize=(10, 10))
    rows = 3
    columns = len(google_scholar_info)
    count  =1
    for idx, img in enumerate(google_scholar_info):
        fig.add_subplot(rows, columns, count)
        count += 1 
        plt.imshow(img)
        plt.axis('off')
    for idx, img in enumerate(wiki_info):
        fig.add_subplot(rows, columns, count)
        count += 1 
        plt.imshow(img)
        plt.axis('off')
    for idx, img in enumerate(open_alex_info):
        fig.add_subplot(rows, columns, count)
        count += 1 
        plt.imshow(img)
        plt.axis('off')
    plt.show()
