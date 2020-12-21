# -*- coding: utf-8 -*-
# Generic front-end code for making biological dataset browsing interfaces.
# Author: Akshay Balsubramani

import dash_core_components as dcc, dash_html_components as html
import app_config
import base64


# ================================================================
# =================== Component styles/layouts ===================
# ================================================================


legend_font_macro = {
    'family': 'sans-serif', 
    'size': app_config.params['legend_font_size'], 
    'color': app_config.params['legend_font_color'] 
}

colorbar_font_macro = {
    'family': 'sans-serif', 
    'size': 8, 
    'color': app_config.params['legend_font_color'] 
}

hm_font_macro = {
    'family': 'sans-serif', 
    'size': 8, 
    'color': app_config.params['legend_font_color'] 
}

style_unselected = {
    'marker': {
        'size': 2.5, 
        'opacity': 1.0
    }
}

style_selected = {
    'marker': {
        'size': 6.0, 
        'opacity': 1.0
    }
}

style_outer_dialog_box = {
    # 'user-select': 'none', '-moz-user-select': 'none', '-webkit-user-select': 'none', '-ms-user-select': 'none', 
    'padding': 10, 
    'margin': 5, 
    # 'borderRadius': 5, 
    'border': 'thin lightgrey solid'
}

style_invis_dialog_box = {
    # 'user-select': 'none', '-moz-user-select': 'none', '-webkit-user-select': 'none', '-ms-user-select': 'none', 
    'padding': 0, 
    'margin': 5
}

style_hm_colorbar = {
    'len': 0.3, 
    'thickness': 20, 
    'xanchor': 'left', 
    'yanchor': 'top', 
    'title': app_config.params['hm_colorvar_name'], 
    'titleside': 'top', 
    'ticks': 'outside', 
    'titlefont': legend_font_macro, 
    'tickfont': legend_font_macro
}

style_text_box = {
    'textAlign': 'center', 
    'width': '100%', 
    'color': app_config.params['font_color']
}

style_legend = {
    'font': legend_font_macro, 
    # bgcolor=app_config.params['legend_bgcolor'], 
    # 'borderwidth': app_config.params['legend_borderwidth'], 
    'padding': 0, 
    'margin': 0, 
    'border': 'thin lightgrey solid', 
    'traceorder': 'normal', 
    'orientation': 'h'
}


def create_hm_layout(
    scatter_frac_domain=0.10, scatter_frac_range=0.08, geneview_range=0.05, 
    show_legend=False, clustersep_coords=None, 
    xaxis_label=True, yaxis_label=True
):
    shape_list = []
    for x in clustersep_coords:
        shape_list.append({
            'type': 'line',
            'x0': x, 'x1': x, 'y0': -0.3, 'y1': 0.3, 'yref': 'y2', 
            'line': { 'color': 'white', 'width': 2 }
        })
    xtext = 'Cell lines' if xaxis_label else ''
    ytext = 'Genes' if yaxis_label else ''
    hm_layout = {
        'annotations': [{
                'x': 0.5, 'y': 1.05, 'showarrow': False, 
                'font': { 'family': 'sans-serif', 'size': 15, 'color': app_config.params['legend_font_color'] }, 
                'text': xtext,
                'xref': 'paper', 'yref': 'paper'
            }, 
            {
                'x': 0.0, 'y': 0.5, 'showarrow': False, 
                'font': { 'family': 'sans-serif', 'size': 15, 'color': app_config.params['legend_font_color'] }, 
                'text': ytext, 
                'textangle': -90, 
                'xref': 'paper', 'yref': 'paper'
            }
        ], 
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 30 }, 
        'clickmode': 'event',  # https://github.com/plotly/plotly.js/pull/2944/
        'hovermode': 'closest', 
        'uirevision': 'Default dataset', 
        'xaxis': {
            'showticklabels': False, 'side': 'top', 
            'tickcolor': app_config.params['legend_bgcolor'], 
            'tickfont': { 'family': 'sans-serif', 'size': app_config.params['hm_font_size'], 'color': app_config.params['legend_font_color'] }, 
            # 'title': {'text': 'Cell lines', 'font': legend_font_macro, 'side': 'top' }, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'domain': [scatter_frac_domain, 1]
        }, 
        'xaxis2': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [0, scatter_frac_domain], 
            'range': [-1, 0.2]
        }, 
        'yaxis': {
            'showticklabels': False, #'side': 'right', 
            'tickcolor': app_config.params['legend_bgcolor'], 
            # 'title': {'text': 'Genes', 'font': legend_font_macro }, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'domain': [0, 1-scatter_frac_range-geneview_range] 
        }, 
        'yaxis2': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [1-scatter_frac_range-geneview_range, 1-geneview_range], 
            'range': [-0.2, 0.5]
        }, 
        'yaxis3': {
            'showticklabels': False, #'side': 'right', 
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [1-geneview_range, 1]
        }, 
        'legend': style_legend, 
        'showlegend': show_legend, 
        'plot_bgcolor': app_config.params['bg_color'], 
        'paper_bgcolor': app_config.params['bg_color'], 
        'shapes': shape_list
    }
    return hm_layout


def create_scatter_layout(annotations):
    return {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 20}, 
        'clickmode': 'event',  # https://github.com/plotly/plotly.js/pull/2944/
        'hovermode': 'closest', 
        'uirevision': 'Default dataset',     # https://github.com/plotly/plotly.js/pull/3236
        'xaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'style': {'display': 'none'}
        }, 
        'yaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'style': {'display': 'none'}
        }, 
        'legend': style_legend, 
        'annotations': annotations, 
        'plot_bgcolor': app_config.params['bg_color'], 
        'paper_bgcolor': app_config.params['bg_color']
    }



# ======================================================
# =================== Component divs ===================
# ======================================================


# Default dataset first in the given list of dataset options.
def create_div_select_dataset(dataset_options):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='three columns', 
                children=[
                    html.A(
                        html.Button(
                            id='download-layout-button', 
                            children='Get TSV', 
                            style=style_text_box, 
                            n_clicks=0, 
                            n_clicks_timestamp=0
                        ), 
                        id='download-layout-link',
                        download="selected_layout.tsv", 
                        href="",
                        target="_blank", 
                        style={
                            'width': '100%', 
                            'textAlign': 'center', 
                            'color': app_config.params['font_color']
                        }
                    )], 
                style={ 'padding-top': '10px' }
            )]
    )


div_landscape_select = html.Div(
    className='row', 
    children=[
        html.Div(
            className='row', 
            children=[
                html.Div(
                    className='six columns', 
                    children=[
                        dcc.Input(
                            id='pointset-name', 
                            type='text', 
                            placeholder="Store gene set with name...", 
                            value='', 
                            style={'width': '100%'}
                        )]
                ), 
                html.Div(
                    className='six columns', 
                    children=[
                        dcc.Dropdown(
                            id='list-pointsets', 
                            placeholder="Select stored gene sets to load", 
                            multi=True
                        )]
                )], style={'display': 'none'}
        )]
)


div_go_panel = html.Div(
    className="row",
    children = [
        dcc.Graph(
            id='goenrich-panel', 
            config={'displaylogo': False, 'displayModeBar': True}
        ), 
        html.Div(
            className='row', 
            children=[
                html.Div(
                    className='nine columns',
                    children=[
                        html.P(
                            "# terms to display: ",
                            style={ 'textAlign': 'right', 'width': '100%', 'color': app_config.params['font_color'] }
                        )], 
                    style={'padding-top': '7px'}
                ),
                html.Div(
                    className='three columns', 
                    children=[
                        dcc.Input(
                            id='select-topk-goterms', 
                            type='text', 
                            value='20', 
                            style={'textAlign': 'center', 'width': '100%'}
                        )]
                )], 
            style={'padding-top': '10px'}
        )]
)


def create_div_hm_panel(point_names):
    return html.Div(
        className="row",
        children = [
            html.Div(
                className='row', 
                children=[
                    html.Div(
                        id='num-selected-counter', 
                        className='four columns', 
                        children='# selected: ', 
                        style={
                            # 'display': 'none', 
                            'textAlign': 'center', 
                            'color': app_config.params['font_color'], 
                            'padding-top': '0px'
                        }
                    ), 
                    html.Div(
                        className='three columns', 
                        children = [
                            dcc.RadioItems(
                                id='select-hm-dataset', 
                                options=[ {'label': v, 'value': v} for v in ['Mutation', 'Expression'] ], 
                                style=legend_font_macro, 
                                labelStyle={
                                    # 'display': 'inline-block', 
                                    'margin-right': '5px'
                                }, 
                                value='Mutation'
                            )]
                    ), 
                    html.Div(
                        className='three columns', 
                        children = [
                            dcc.Dropdown(
                                id='select-geneview', 
                                placeholder="View gene", #multi=True, 
                                style={'display': 'inline-block', 'width': '100%', 'textAlign': 'center'}
                            )]
                    )
                ], 
                style={'margin': 5}
            ), 
            dcc.Loading(
                id="loading-heatmap", 
                children=[
                    dcc.Graph(
                        id='main-heatmap', 
                        style={ 'height': '60vh' }, 
                        config={'displaylogo': False, 'displayModeBar': True}
                    )], type="cube"
            )]
    )

# ==================================================================
# =================== Aggregating component divs ===================
# ==================================================================

import numpy as np
# go_termIDs = np.load(app_config.params['gotermIDs_path'])
# go_termnames = np.load(app_config.params['gotermnames_path'])

def create_div_mainctrl(
    point_names, feat_names, more_colorvars, cancer_types, upload_asset, download_asset
):
#     download_image = app_config.params['download_img_path']#,
#     encoded_image = base64.b64encode(open(download_image, 'rb').read())
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='row', 
                children=[
                    html.Div(
                        className='four columns', 
                        children=[
                            dcc.Dropdown(
                                id='gene-selection', 
                                options = [ {'value': point_names[i], 'label': "{}".format(point_names[i])} for i in range(len(point_names)) ], 
                                placeholder="Gene...", multi=True, 
                                style={'height': '55px', 'display': 'inline-block', 'width': '100%', 'textAlign': 'center'}
                            )], 
                        style={'fontSize': 12, 'margin': 5}
                    ), 
                    html.Div(
                        className='four columns', 
                        id='div-cell-line-lookup', 
                        children=[
                            dcc.Dropdown(
                                id='cell-line-lookup', 
                                options = [{'value': n, 'label': n} for n in more_colorvars] + [{'value': gn, 'label': ' '.join(gn.split('_'))} for gn in feat_names], 
                                value=app_config.params['default_color_var'], 
                                placeholder="Cell line...", 
                                clearable=True, 
                                style={'white-space':'nowrap', 'text-overflow': 'ellipsis', 
                                       'height': '55px', 'display': 'inline-block', 'width': '100%', 'textAlign': 'center' }
                            )], 
                        style={'fontSize': 12, 'margin': 5}, 
                        n_clicks_timestamp=0
                    ), 
                    html.Div(
                        className='four columns', 
                        id='div-tissue-lookup', 
                        children=[
                            dcc.Dropdown(
                                id='tissue-type-lookup', 
                                options = [{'value': n, 'label': n} for n in cancer_types], 
                                placeholder="Tissue type...", 
                                value='', 
                                style={'height': '55px', 'display': 'inline-block', 'width': '100%', 'textAlign': 'center'}
                            )], 
                        style={'fontSize': 12, 'margin': 5}, 
                        n_clicks_timestamp=0
                    )], 
                style={}
            ), 
            html.Div(
                className='row', 
                children=[
                    html.Div(
                        className='five columns', 
                        children=[
                            html.Div(
                                className='row', 
                                children='Gene set download/upload:', 
                                style={
                                    'textAlign': 'center', 
                                    'color': app_config.params['font_color'], 
                                    'padding-top': '0px'
                                }
                            ), 
                            html.Div(
                                className='row', 
                                children=[
                                    html.Div(
                                        className='six columns', 
                                        children=[
                                            html.A(
                                                html.Button(
                                                    id='download-button', 
                                                    children='Save', 
                                                    style=style_text_box, 
                                                    n_clicks=0, 
                                                    n_clicks_timestamp=0
                                                ), 
                                                id='download-set-link',
                                                download="selected_set.csv", 
                                                href="",
                                                target="_blank", 
                                                style={
                                                    'width': '100%', 
                                                    'textAlign': 'center', 
                                                    'color': app_config.params['font_color']
                                                }
                                            )], 
                                        style={'padding-top': '0px'}
                                    ), 
                                    html.Div(
                                        className='six columns', 
                                        children=[
                                            dcc.Upload(
                                                id='upload-pointsets',
                                                children=html.Div([
                                                    html.Button(
                                                        id='upload-button', 
                                                        children='Load', 
                                                        style=style_text_box, 
                                                        n_clicks=0, 
                                                        n_clicks_timestamp=0
                                                    ) # html.Img( src=upload_asset, #'data:image/png;base64,{}'.format(encoded_image), style={ 'height' : '40%', 'width' : '40%', 'float' : 'right', 'position' : 'relative', 'padding-top' : 0, 'padding-right' : 0 })
                                                ]),
                                                style={
                                                    'width': '100%', 
                                                    'textAlign': 'center', 
                                                    'color': app_config.params['font_color']
                                                }, 
                                                multiple=True
                                            )]
                                    )]
                            )], style={ 'border': 'thin lightgrey solid',  'margin': 5 }
                    ), 
                    html.Div(
                        # id='num-selected-counter', 
                        className='five columns', 
                        children=[
                            html.Button(
                                id='hm-button', 
                                children='Plot essentiality >>', 
                                style=style_text_box, 
                                n_clicks=0, 
                                n_clicks_timestamp=0
                            )], 
                        style={'padding-top': '0px'}
                    )
                ], style={}
            )
        ]
    )



"""
Main layout.
"""
import numpy as np

def create_div_mainapp(point_names, feat_names, cancer_types, upload_asset, download_asset, more_colorvars=[]):
    browser_div = html.Div(
        className="browser-div", 
        children=[
            html.Div(
                className="row", 
                children=[
                    html.Div(
                        className="seven columns",
                        children=[
                            create_div_mainctrl(point_names, feat_names, more_colorvars, cancer_types, upload_asset, download_asset), 
                            dcc.Graph(
                                id='landscape-plot',
                                config={'displaylogo': False, 'displayModeBar': True}, 
                                style={ 'height': '100vh'}
                            ), 
                            create_div_select_dataset(app_config.params['dataset_options']), 
                            div_landscape_select
                        ]
                    ), 
                    html.Div(
                        className='five columns', 
                        children=[
                            create_div_hm_panel(point_names), 
                            div_go_panel
                        ], 
                        style=style_invis_dialog_box
                    )]
            ), 
            html.Div([ html.Pre(id='test-select-data', style={ 'color': app_config.params['font_color'], 'overflowX': 'scroll' } ) ]),     # For testing purposes only!
            html.Div(
                className='row', 
                children=[ 
                    dcc.Markdown(
                        """Source [repository](https://github.com/kundajelab/coessentiality-browser). """
                         + """Queries? [Contact](abalsubr@stanford.edu). """
                    )], 
                style={
                    'textAlign': 'center', 
                    'color': app_config.params['font_color'], 
                    'padding-bottom': '10px'
                }
            )],
        style={ 
            'width': '100vw', 
            'max-width': 'none'
        }
    )
    return html.Div(
        className="container", 
        children=[
            html.Div(
                className='row', 
                children=[
                    html.H1(
                        id='title', 
                        children=app_config.params['title'], 
                        style=style_text_box
                    )]
            ), 
            browser_div
        ],
        style={ 
            'width': '100vw', 
            'max-width': 'none'
        }
    )
    
