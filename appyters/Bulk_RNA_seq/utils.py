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
import warnings
import numpy as np

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
        
        


robjects.r('''limma <- function(rawcount_dataframe, design_dataframe, filter_genes=FALSE, adjust="BH") {
    #print(rawcount_dataframe)
    #print(design_dataframe)
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
    return(limma_dataframe)
}
''')