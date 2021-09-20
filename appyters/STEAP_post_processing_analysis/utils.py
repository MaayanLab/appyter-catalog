from IPython.display import HTML, display, Markdown, IFrame, FileLink
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models import ColumnDataSource
from bokeh.plotting import show
import pandas as pd
import requests
from bs4 import BeautifulSoup
import re
from zipfile import ZipFile


def display_bokeh_df(bokeh_df, counter, caption):
    bokeh_columns = [TableColumn(field=c, title=c) for c in bokeh_df.columns] # bokeh columns
    bokeh_datatable = DataTable(
        columns=bokeh_columns,
        source=ColumnDataSource(bokeh_df),
        width=950
    ) # bokeh table
    show(bokeh_datatable)
    counter = make_caption(counter, caption, 'Table')
    return counter

def make_caption(counter, caption, image):
    display(Markdown(f"*{image} {counter}. {caption}*"))
    counter += 1
    return counter

def display_link(url):
    return HTML(f"<a href=\"./{url}\" target='_blank'>Download {url}</a>")
    
def create_download_link(df, filename):
    df.to_excel(filename, index=False)
    return HTML(f"<a href=\"./{filename}\" target='_blank'>Download {filename}</a>")

def create_zip_file(outfile, infiles):
    with ZipFile(outfile, 'w') as zipObj:
        for f in infiles:
            zipObj.write(f)

def get_esmu_data():
    url = 'https://github.com/erwinerdem/STEAP/tree/master/esmu/'
    r = requests.get(url)
    soup = BeautifulSoup(r.text, features="lxml")
    scrna_data = soup.body.findAll(text=re.compile('.csv'))
    scrna_dict = {
        d.split('.')[0]: pd.read_csv("https://github.com/erwinerdem/STEAP/raw/master/esmu/"+d ,index_col=0)
        for d in scrna_data
    }
    return scrna_dict