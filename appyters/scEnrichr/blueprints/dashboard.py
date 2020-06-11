def dashboard(app, url_prefix='/dashboard', DATA_DIR=''):
  import os
  import json
  import pandas as pd
  from flask import Blueprint
  from functools import lru_cache
  import uuid
  import dash
  import dash_auth
  import dash_table as dt
  import dash_core_components as dcc
  import dash_html_components as html
  from dash.dependencies import Input, Output
  from dash.exceptions import PreventUpdate
  from werkzeug.utils import secure_filename

  dashapp = dash.Dash(
    'dashboard',
    server=app,
    routes_pathname_prefix=url_prefix + '/',
  )

  def sanitize_uuid_path(val):
    return str(uuid.UUID(val[len(url_prefix) + 1:]))

  @lru_cache()
  def df_for_session(session):
    return pd.read_csv(
      os.path.join(DATA_DIR, session, 'df.tsv'),
      sep='\t'
    )

  @lru_cache()
  def df_umap_for_session(session):
    return pd.read_csv(
      os.path.join(DATA_DIR, session, 'df_umap.tsv'),
      sep='\t'
    )

  @lru_cache()
  def df_enrich_for_session(session):
    return pd.read_csv(
      os.path.join(DATA_DIR, session, 'df_enrich.tsv'),
      sep='\t'
    )

  locks = {}

  def lock_for_session(session):
    if session not in locks:
      locks[session] = False
    return locks[session]

  prevClickDatas = {}
  def prevClickData_for_session(session):
    if session not in prevClickDatas:
      prevClickDatas[session] = None
    return prevClickDatas[session]

  def umap_figure_for_session(session):
    df_umap = df_umap_for_session(session)
    return {
      'data': [
        dict(
            x=df_umap[df_umap['Cluster'] == cluster]['UMAP-1'],
            y=df_umap[df_umap['Cluster'] == cluster]['UMAP-2'],
            mode='markers',
            opacity=1,
            marker=dict(
                size=8,
                line=dict(width=0.5, color='white'),
            ),
            text=df_umap[df_umap['Cluster'] == cluster]['Cluster'],
            name='Cluster {} '.format(cluster),
        ) for cluster in sorted(df_umap['Cluster'].unique())
      ],
      'layout': dict(
          xaxis=dict(
              title='UMAP-1',
          ),
          yaxis=dict(
              title='UMAP-2',
          ),
          margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
          legend={'x': 0, 'y': 1},
          hovermode='closest'
      ),
    }

  dashapp.layout = html.Div(children=[
      dcc.Location(id='url', refresh=False),
      dcc.Graph(id='umap'),
      html.H2(id='cluster-header'),
      html.Label(id='enrichr-link'),
      dt.DataTable(
          id='data-table',
          columns=[
              {'name': 'cluster', 'id': 'cluster'},
              {'name': 'rank', 'id': 'rank'},
              {'name': 'direction', 'id': 'direction'},
              {'name': 'term', 'id': 'term'},
              {'name': 'category', 'id': 'category'},
              {'name': 'pvalue', 'id': 'pvalue', 'type': 'numeric', 'format': { 'specifier': '.3' } },
              {'name': 'library', 'id': 'library'},
          ],
          sort_action='native',
          sort_mode='multi',
          filter_action='native',
          filter_query='{pvalue} < 0.05 && {direction} = up',
          page_action='native',
          sort_by=[{ 'column_id': 'pvalue', 'direction': 'asc' }],
          style_as_list_view=True,
          style_header={
              'backgroundColor': 'rgb(200, 200, 200)',
              'fontWeight': 'bold'
          },
          style_table={
              'overflowX': 'auto',
              'width': '100%',
              'minWidth': '100%',
          },
          style_data={
              'whiteSpace': 'normal',
              'height': 'auto',
          },
          css=[
              {
                  'selector': '.dash-cell div.dash-cell-value',
                  'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;',
              },
          ],
          style_data_conditional=[
              {
                  'if': {'row_index': 'odd'},
                  'backgroundColor': 'rgb(230, 230, 230)'
              }
          ],
          style_cell_conditional=[
              {
                  'if': { 'column_id': 'term' },
                  'textAlign': 'left',
                  'maxWidth': '20vw',
              },
              {
                  'if': { 'column_id': 'category' },
                  'textAlign': 'left',
              },
              {
                  'if': { 'column_id': 'direction' },
                  'textAlign': 'center',
              },
              {
                  'if': { 'column_id': 'library' },
                  'textAlign': 'left',
                  'maxWidth': '10vw',
              },
          ],
      ),
  ])

  @dashapp.callback(dash.dependencies.Output('umap', 'figure'),
                [dash.dependencies.Input('url', 'pathname')])
  def load_data(pathname):
    session = sanitize_uuid_path(pathname)
    return umap_figure_for_session(session)


  lock = False
  prevClickData = None

  @dashapp.callback(
      [Output('cluster-header', 'children'), Output('enrichr-link', 'children'), Output('data-table', 'data')],
      [dash.dependencies.Input('url', 'pathname'), Input('umap', 'clickData'), Input('umap', 'hoverData')]
  )
  def update_click(pathname, clickData, hoverData):
    session = sanitize_uuid_path(pathname)
    df_umap = df_umap_for_session(session)
    df_enrich = df_enrich_for_session(session)
    lock = lock_for_session(session)
    prevClickData = prevClickData_for_session(session)
    # Initial state
    if not clickData and not hoverData:
        return ['Click to cluster to select', '', []]
    # Get relevant evt
    if prevClickData != clickData: # Click
        lock = not lock
        prevClickData = clickData
        evt = clickData
    elif lock: # Hover but locked
        raise PreventUpdate
    else: # Hover and not locked
        evt = hoverData
    # Get cluster
    cluster = evt['points'][0]['text']
    matches = df_enrich[df_enrich['cluster'] == cluster]
    if matches.size == 0:
        return [
            'Cluster {} ({} samples)'.format(cluster, df_umap[df_umap['Cluster'] == cluster].shape[0]),
            'No data for this cluster',
            [],
        ]
    # Update
    link = matches.iloc[0]['link']
    data = matches.to_dict('records')
    return [
        'Cluster {} ({} samples)'.format(cluster, df_umap[df_umap['Cluster'] == cluster].shape[0]),
        ['Enrichr Link for Cluster ', html.A(link, href=link)],
        data,
    ]
