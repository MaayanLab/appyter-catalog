#import packages
from time import sleep
from collections import defaultdict
from IPython.display import display,FileLink, HTML, Markdown, IFrame
import requests
import plotly.express as px
import plotly
from scholarly import scholarly
import datetime
import math
from bs4 import BeautifulSoup
import plotly.io as pio
pio.renderers.default = 'notebook'

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
    display(Markdown("*Figure {}. {}*".format(counter, caption)))
    if title != None:
        display(FileLink("output_images/"+title+'.png' , result_html_prefix=str('Download Figure {} (PNG): '.format(counter))))
    counter += 1
    return counter


def display_summary_text_from_google_scholar(affiliation_from_google_scholar = None, h_index_from_google_scholar = None, interests_from_google_scholar = [], h_index_from_google_scholar_last_5 = None, ar_index = None, total_times_cited = None, name_of_researcher = None):

    # <h3 style="color: #333;">Organization Affiliation from RePORTER</h3>
    # <p style="color: #333;">{organization_from_reporter}</p>
    display(Markdown("## Text Summary Information for {} (From Google Scholar) ##".format(name_of_researcher)))
    html_for_text_display = ""
    html_for_text_display += "<div style=\"display: grid; grid-template-columns: 200px 200px 200px; grid-gap: 20px; padding: 20px; width: 750px\">\n"
    if affiliation_from_google_scholar != None:
        html_for_text_display += f"""
    <div style="background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;">
        <h3 style="color: #333;">Institution Affiliation</h3>
        <p style="color: #333;">{affiliation_from_google_scholar}</p>
    </div>
    """
    if h_index_from_google_scholar != None:
        html_for_text_display += f"""
        <div style="background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;">
        <h3 style="color: #333;">H-index: {h_index_from_google_scholar}</h3>
        <h3 style="color: #333;">H-index over last 5 years: {h_index_from_google_scholar_last_5} </h3>
        <h3 style="color: #333;">Age Related Index: {int(ar_index)}</h3>
        <h3 style="color: #333;">Total Times Cited: {total_times_cited}</h3>
        
    </div>
    """
    if len(interests_from_google_scholar) > 0:
        html_for_text_display += "<div style=\"background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;\">\n"
        html_for_text_display += "<h3 style=\"color: #333;\">Research Interests</h3>\n"
        html_for_text_display += "<ul>\n"
        for interest in interests_from_google_scholar:
            html_for_text_display += f"<li><p style=\"color: #333;\">{interest}</p> </li>\n"
        html_for_text_display += "</ul>\n"
        html_for_text_display += "</div\n"
    html_for_text_display += "</div>"
    display(HTML(html_for_text_display))

def query_google_citation(name_of_researcher):
    search_query = scholarly.search_author(name_of_researcher)
    try:
        author  = next(search_query)
        display(Markdown("### Link to [Google Scholar Page](https://scholar.google.com/citations?user={}) for {}".format(author['scholar_id'], name_of_researcher)))
        summary_info = scholarly.fill(author, sections=['counts', 'indices'])
        # print(summary_info)
        # Print the titles of the author's publications
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
        display_summary_text_from_google_scholar(affiliation_from_google_scholar, h_index_from_google_scholar, interests_from_google_scholar, h_index_from_google_scholar_last_5, ar_index, total_times_cited, name_of_researcher)
        return citation_dict
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
            # print(id_of_researcher)
            # print(final_name)
            url_for_papers_final = f"https://api.semanticscholar.org/graph/v1/author/{id_of_researcher}?fields=name,citationCount,paperCount,hIndex,aliases,papers.year,papers.citationCount"
            citation_response = requests.get(url_for_papers_final)
            if citation_response.status_code == 200:
                display(Markdown("### Link to [Semantic Scholar Page](https://www.semanticscholar.org/author/{}/{}) for {}".format(final_name, id_of_researcher, name_of_researcher)))
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

def display_summary_text_from_wikipedia(wiki_institution = '', wiki_known_for = [], wiki_field_interests = [], wiki_awards = [], name_of_researcher = None):

    display(Markdown("## Text Summary Information for {} (From Wikipedia) ##".format(name_of_researcher)))
    html_for_text_display = ""
    html_for_text_display += "<div style=\"display: grid; grid-template-columns: 150px 150px 150px 150px; grid-gap: 20px; padding: 20px; width: 750px\">\n"
    if wiki_institution != '':
        html_for_text_display += f"""
    <div style="background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 150px;">
        <h3 style="color: #333;">Institution Affiliation</h3>
        <p style="color: #333;">{wiki_institution}</p>
    </div>
    """
    if len(wiki_awards) > 0:
        html_for_text_display += "<div style=\"background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 150px;\">\n"
        html_for_text_display += "<h3 style=\"color: #333;\">Awards</h3>\n"
        html_for_text_display += "<ul>\n"
        for award in wiki_awards:
            html_for_text_display += f"<li><p style=\"color: #333;\">{award}</p> </li>\n"
        html_for_text_display += "</ul>\n"
        html_for_text_display += "</div>\n"

    if len(wiki_known_for) > 0:
        html_for_text_display += "<div style=\"background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 150px;\">\n"
        html_for_text_display += "<h3 style=\"color: #333;\">Known For</h3>\n"
        html_for_text_display += "<ul>\n"
        for known in wiki_known_for:
            html_for_text_display += f"<li><p style=\"color: #333;\">{known}</p> </li>\n"
        html_for_text_display += "</ul>\n"
        html_for_text_display += "</div>\n"

    if len(wiki_field_interests) > 0:
        html_for_text_display += "<div style=\"background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 150px;\">\n"
        html_for_text_display += "<h3 style=\"color: #333;\">Research Fields</h3>\n"
        html_for_text_display += "<ul>\n"
        for interest in wiki_field_interests:
            html_for_text_display += f"<li><p style=\"color: #333;\">{interest}</p> </li>\n"
        html_for_text_display += "</ul>\n"
        html_for_text_display += "</div>\n"
    html_for_text_display += "</div>"
    display(HTML(html_for_text_display))

def display_summary_text_from_openalex(institution = '', interests = [], h_index = None, i10_index = None, total_times_cited = None, name_of_researcher = None):

    display(Markdown("## Text Summary Information for {} (From OpenAlex API) ##".format(name_of_researcher)))
    html_for_text_display = ""
    html_for_text_display += "<div style=\"display: grid; grid-template-columns: 200px 200px 200px; grid-gap: 20px; padding: 20px; width: 750px\">\n"
    if institution != '':
        html_for_text_display += f"""
    <div style="background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;">
        <h3 style="color: #333;">Institution Affiliation</h3>
        <p style="color: #333;">{institution}</p>
    </div>
    """

    if h_index != None:
        html_for_text_display += f"""
        <div style="background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;">
        <h3 style="color: #333;">H-index: {h_index}</h3>
        <h3 style="color: #333;">I10-index : {i10_index} </h3>
        <h3 style="color: #333;">Total Times Cited: {total_times_cited}</h3>
        
    </div>
    """
    if len(interests) > 0:
        html_for_text_display += "<div style=\"background-color: #f1f1f1; padding: 5px; border-radius: 5px; width: 200px;\">\n"
        html_for_text_display += "<h3 style=\"color: #333;\">Awards</h3>\n"
        html_for_text_display += "<ul>\n"
        for interest in interests:
            html_for_text_display += f"<li><p style=\"color: #333;\">{interest}</p> </li>\n"
        html_for_text_display += "</ul>\n"
        html_for_text_display += "</div>\n"

    html_for_text_display += "</div>"
    display(HTML(html_for_text_display))


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
                for concept in first_result['x_concepts'][:6]:
                    interests.append(concept['display_name'])
            
            display_summary_text_from_openalex(institution,interests, h_index,i10_index, total_times_cited, name_of_researcher)
        else:
            print("No matching information for this researcher from OpenAlex")
    else:
        print("Error in querying from OpenAlex API.")
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
                display(Markdown("### Link to [Wiki Page](https://en.wikipedia.org/wiki/{}) for {}".format(title_to_get_content_from.replace(" ","_"), name_of_researcher)))
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
                                            if index != len(td.contents):
                                                wiki_institution += content.text + " and "
                                            else:
                                                wiki_institution += content.text
                                elif "Awards" in th.text:
                                    for content in td.contents:
                                        if type(content) == type(BeautifulSoup("", "html.parser").new_tag("tag")) and len(content.text) > 0:
                                            wiki_awards.append(content.text)
                                            # print(content.text)
                        display_summary_text_from_wikipedia(wiki_institution, wiki_known_for, wiki_field_interests, wiki_awards, name_of_researcher)
                    else:
                        print("No data to show")
                except:
                    print("Error in reading the content from Wikipedia page")
            else:
                print("No matching Wikipedia Page for this researcher")
        else:
            print("No matching Wikipedia Page for this researcher")
    else:
        print("Error in querying Wikipedia for the name of researcher")         
