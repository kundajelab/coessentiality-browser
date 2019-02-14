# -*- coding: utf-8 -*-
# Generic front-end code for making biological dataset browsing interfaces.
# Author: Akshay Balsubramani

import dash_core_components as dcc, dash_html_components as html
import app_config


# ========================================================
# =================== Component styles ===================
# ========================================================


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
    'padding': 5, 
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

style_upload = {
    'width': '100%', 
    'border': 'thin lightgrey solid',
    'textAlign': 'center', 
    'color': app_config.params['font_color'], 
    'padding-top': '5px', 
    'padding-bottom': '5px'
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


def create_hm_layout(scatter_frac_domain):
    hm_layout = {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 0 }, 
        'clickmode': 'event+select',  # https://github.com/plotly/plotly.js/pull/2944/
        'hovermode': 'closest', 
        'uirevision': 'Default dataset', 
        'xaxis': {
            'automargin': True, 
            'showticklabels': True, 
            'side': 'top', 
            'tickcolor': app_config.params['legend_bgcolor'], 
            # 'tickangle': 60, 
            'tickfont': {
                'family': 'sans-serif', 
                'size': app_config.params['hm_font_size'], 
                'color': app_config.params['legend_font_color'] 
            }, 
            'dtick': 1, 
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [scatter_frac_domain, 1]
        }, 
        'xaxis2': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [0, scatter_frac_domain]
        }, 
        'yaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False 
        }, 
        'legend': style_legend, 
        'plot_bgcolor': app_config.params['bg_color'], 
        'paper_bgcolor': app_config.params['bg_color']
    }
    return hm_layout


def create_scatter_layout(annotations):
    return {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 0}, 
        'clickmode': 'event+select',  # https://github.com/plotly/plotly.js/pull/2944/
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


def create_div_plotcolor(feat_names, more_colorvars):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[ html.P('Plot as color:') ], 
                style={
                    'textAlign': 'center', 
                    'color': app_config.params['font_color'], 
                    'padding-top': '5px'
                }
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.Dropdown(
                        id='landscape_color', 
                        options = [{
                            'value': app_config.params['default_color_var'], 
                            'label': app_config.params['default_color_var']
                        }] + [
                            {'value': n, 'label': n} for n in more_colorvars
                        ] + [
                            {'value': gn, 'label': gn} for gn in feat_names 
                        ], 
                        value=app_config.params['default_color_var'], 
                        placeholder="Select colors to plot", 
                        clearable=False
                    )], 
                style={'padding-top': '0px'}
            )], 
        style=style_invis_dialog_box
    )


def create_div_annot(point_names):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[ html.P('Look up gene(s) in plot:') ], 
                style={
                    'textAlign': 'center', 
                    'color': app_config.params['font_color'], 
                    'padding-top': '5px'
                }
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.Dropdown(
                        id='points_annot', 
                        options = [ {'value': gn, 'label': gn} for gn in point_names ], 
                        placeholder="Select genes", multi=True
                    )], 
                style={'padding-top': '0px'}
            )]
    )


def create_div_align_selection(options_list):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[ 
                    html.Button(
                        id='align-button', 
                        children='Display alignment', 
                        n_clicks=0, 
                        n_clicks_timestamp=0,
                        style=style_text_box
                    )], 
                style={'padding-top': '5px'}
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.RadioItems(
                        id='align-method-selection', 
                        options=[ {'label': v, 'value': v} for v in options_list ], 
                        style=legend_font_macro, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }, 
                        value='Unaligned'
                    )], 
                style={'padding-top': '10px'}
            )]
    )


# Default dataset first in the given list of dataset options.
def create_div_select_dataset(dataset_options):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[
                    html.P(
                        "Browse layout: ", 
                        style=style_text_box
                    )], 
                style={'padding-top': '10px'}
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.Dropdown(
                        id='sourcedata-select', 
                        options = [ {'value': dn, 'label': dn} for dn in dataset_options ], # style={'height': '30px'}, 
                        value=dataset_options[0]
                    )]
            )], 
        style=style_outer_dialog_box
    )


# Synchronizing landscape plot with heatmap.
def create_div_sync_selection():
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='five columns', 
                children=[
                    html.Button(
                        id='hm-highlight-button', 
                        children='<< Highlight from heatmap', 
                        style=style_text_box, 
                        n_clicks='0', 
                        n_clicks_timestamp='0'
                    )], style={'padding-top': '5px'}
            ), 
            html.Div(
                className='three columns', 
                children=[
                    dcc.Checklist(
                        id='toggle-hm-zoom', 
                        options=[
                            {'label': 'Zoom', 'value': 'on'}
                        ],
                        values=[], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }
                    )], 
                style={'padding-top': '0px'}
            ), 
            html.Div(
                className='four columns', 
                children=[
                    dcc.RadioItems(
                        id='main-heatmap-roworder', 
                        options=[ 
                            {'label': 'Cocluster', 'value': 'Cocluster'}, 
                            {'label': 'Sort by color', 'value': 'Sort by color'}], 
                        style=legend_font_macro, 
                        labelStyle={
                            # 'display': 'inline-block', 
                            'margin-right': '5px'
                        }, 
                        value='Sort by color'
                    )]
            )
        ], 
        style=style_outer_dialog_box
    )




# ================================================
# ==================== Layout ====================
# ================================================


def create_div_go_ctrl(point_names):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='six columns', 
                children=[
                    dcc.Dropdown(
                        id='geneset-select', 
                        options = [ {'value': gn, 'label': gn} for gn in point_names ], # style={'height': '30px'}, 
                        value = [], 
                        placeholder="Gene set for enrichment", multi=True
                    )], 
                style={'padding-top': '0px'}
            ), 
            html.Div(
                className='three columns',
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
                    )
                ]
            )], 
        style=style_outer_dialog_box
    )


div_ctrl_landscape = html.Div(
    className='row', 
    children=[
        html.Div(
            className='four columns', 
            children=[ 
                html.P(
                    'Select genes from:', 
                    style=style_text_box
                )], 
            style={'padding-top': '5px'}
        ),  
        html.Div(
            className='four columns', 
            children=[
                dcc.RadioItems(
                    id='selection-mode-radio', 
                    options=[
                        {'label': 'Scatterplot', 'value': 'scatter'}, 
                        {'label': 'Heatmap', 'value': 'heatmap'}
                    ], 
                    style=legend_font_macro, 
                    labelStyle={
                        'display': 'inline-block', 
                        'margin-right': '5px'
                    }, 
                    value='scatter'
                )], 
            style={'padding-top': '5px'}
        ), 
        html.Div(
            className='four columns', 
            children=[ 
                html.Button(
                    id='plot-selection-button-old', 
                    children='old button', 
                    n_clicks=0, 
                    n_clicks_timestamp=0,
                    style=style_text_box
                )], 
            style={'padding-top': '5px'}
        )]
)


div_select_points = html.Div(
    className='row', 
    children=[
        html.Div(
            id='num-selected-counter', 
            className='three columns', 
            children='# selected: ', 
            style={
                'textAlign': 'center', 
                'color': app_config.params['font_color'], 
                'padding-top': '10px'
            }
        ), 
        html.Div(
            className='five columns', 
            children=[
                dcc.Input(
                    id='pointset-name', 
                    type='text', 
                    placeholder="Store gene set with name...", 
                    value='', 
                    style={'width': '100%'}
                )], 
            style={'padding-top': '5px'}
        ), 
        html.Div(
            className='four columns', 
            children=[
                dcc.Checklist(
                    id='store-status', 
                    options=[
                        {'label': 'Store subset', 'value': 'store'}
                    ],
                    values=[], 
                    style={
                        'textAlign': 'center', 
                        'width': '80%', 
                        'color': app_config.params['font_color']
                    }
                )], 
            style={'padding-top': '10px'}
        )
    ]
)


div_serialization = html.Div(
    className='row', 
    children=[
        html.Div(
            className='four columns', 
            children=[
                html.A(
                    'Save gene set selection',
                    id='download-set-link',
                    download="selected_set.json", 
                    href="",
                    target="_blank", 
                    style={
                        'textAlign': 'center', 
                        'width': '100%', 
                        #'padding-top': '0px', 
                        'color': app_config.params['font_color']
                    }
                )], 
        ), 
        html.Div(
            className='eight columns', 
            children=[
                dcc.Upload(
                    id='upload-pointsets',
                    children=html.Div([
                        'Load gene set from file(s)...'
                    ]),
                    style=style_upload, 
                    multiple=True
                )], 
            style={'padding-top': '0px'}
        )], 
    style=style_invis_dialog_box
)


div_list_pointsets = html.Div(
    className='row', 
    children=[
        html.Div(
            className='eight columns', 
            children=[
                dcc.Dropdown(
                    id='list-pointsets', 
                    placeholder="Select stored gene sets to load", 
                    multi=True
                )], 
            style={'padding-top': '5px'}
        ), 
        html.Div(
            className='four columns', 
            children=[
                dcc.Checklist(
                    id='load-status', 
                    options=[
                        {'label': 'Load subset', 'value': 'load'}
                    ],
                    values=[], 
                    style={
                        'textAlign': 'center', 
                        'width': '80%', 
                        'color': app_config.params['font_color']
                    }
                )], 
            style={'padding-top': '10px'}
        )]
)


div_reviz_scatter = html.Div(
    className='row', 
    children=[
        html.Div(
            className='four columns', 
            children=[ 
                dcc.Checklist(
                    id='reviz-status', 
                    options=[
                        {'label': 'Visualize selection', 'value': 'viz'}
                    ],
                    values=[], 
                    style={
                        'textAlign': 'center', 
                        # 'width': '80%', 
                        'color': app_config.params['font_color']
                    }
                )], 
            style={'padding-top': '10px'}
        ), 
        html.Div(
            className='eight columns', 
            children=[
                dcc.RadioItems(
                    id='reviz-method-selection', 
                    options=[ {'label': v, 'value': v} for v in ['dummy', 'UMAP'] ], 
                    style=legend_font_macro, 
                    labelStyle={
                        'display': 'inline-block', 
                        'margin-right': '5px'
                    }, 
                    value='dummy'
                )], 
            style={'padding-top': '10px'}
        )]
)


div_landscape_select = html.Div(
    className='row', 
    children=[
        div_select_points, 
        div_list_pointsets, 
        div_serialization], 
    style=style_outer_dialog_box
)


div_go_panel = html.Div(
    className="row",
    children = [
        dcc.Graph(
            id='goenrich-panel', 
            # style={ 'height': '30vh' }, 
            config={'displaylogo': False}
        )]
)


def create_div_landscape_ctrl():
    return html.Div(
        className='row', 
            children=[
            dcc.Checklist(
                id='toggle-landscape-whiten', 
                options=[
                    {'label': 'Whiten', 'value': 'on'}
                ],
                values=[], 
                style={
                    'textAlign': 'left', 
                    'width': '80%', 
                    'color': app_config.params['font_color']
                }
            )], 
        style={'padding-top': '0px'}
    )


def create_div_hm_panel():
    return html.Div(
        className="row",
        children = [
            html.Div(
                id='hm-feat-control', 
                children=[]
            ), 
            dcc.Graph(
                id='main-heatmap', 
                style={ 'height': '55vh' }, 
                config={'displaylogo': False}
            )]
    )


# Default dataset first in the given list of dataset options.
def create_div_cosmetic_panel():
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='two columns', 
                children=[
                    html.P(
                        "Interface options: ", 
                        style=style_text_box
                    )], 
                style={'padding-top': '10px'}
            ), 
            html.Div(
                className='three columns', 
                children=[
                    dcc.Slider(
                        id='slider-bg-marker-size-factor', # marks={ 4: {'label': '4'} }
                        min=0, max=8, step=0.2, 
                        value=app_config.params['bg_marker_size_factor']
                    ), 
                    html.Div(
                        id='display-bg-marker-size-factor', 
                        style={
                            'textAlign': 'center', 
                            'color': app_config.params['font_color'], 
                            'padding-top': '10px'
                        }
                    )]
            ), 
            html.Div(
                className='three columns', 
                children=[
                    dcc.Slider(
                        id='slider-marker-size-factor',
                        min=0, max=8, step=0.2, 
                        value=app_config.params['marker_size_factor']
                    ), 
                    html.Div(
                        id='display-marker-size-factor', 
                        style={
                            'textAlign': 'center', 
                            'color': app_config.params['font_color'], 
                            'padding-top': '10px'
                        }
                    )]
            ), 
            html.Div(
                className='four columns', 
                children=[
                    dcc.Checklist(
                        id='toggle-hm-feat-panels', 
                        options=[
                            {'label': 'Per-cluster', 'value': 'bars'}, 
                            {'label': 'Dendrogram', 'value': 'dendrogram'}, 
                            {'label': 'Heatmap', 'value': 'heatmap'}
                        ],
                        values=[], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }
                    )], 
                style={'padding-top': '0px'}
            ), 
            html.Div(
                className='row', 
                children=[
                    dcc.Checklist(
                        id='toggle-debug-panels', 
                        options=[
                            {'label': 'Debug panel', 'value': 'debug-panel'}
                        ],
                        values=[], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }
                    )], 
                style={'padding-top': '0px'}
            )
        ], 
        style=style_outer_dialog_box
    )


def create_div_landscapes(point_names, feat_names, more_colorvars):
    return html.Div(
        className="seven columns",
        children=[
            # div_ctrl_landscape,
            dcc.Graph(
                id='landscape-plot',
                config={'displaylogo': False, 'displayModeBar': True},
		style={ 'height': '100vh'}
            ), 
            create_div_select_dataset(app_config.params['dataset_options']), 
            div_landscape_select
        ]
    )


def create_div_sidepanels(point_names, feat_names, more_colorvars):
    return html.Div(
        className='five columns', 
        children=[
            create_div_landscape_ctrl(), 
            create_div_sync_selection(), 
            create_div_hm_panel(), 
            create_div_plotcolor(feat_names, more_colorvars),
            create_div_annot(point_names), 
            create_div_go_ctrl(point_names), 
            div_go_panel
            # div_reviz_scatter
        ],
        style=style_invis_dialog_box
    )



"""
Main layout.
"""

def create_div_mainapp(point_names, feat_names, more_colorvars=[]):
    return html.Div(
        className="container", 
        children=[
            html.Div(
                className='ten columns offset-by-one', 
                children=[
                    html.H1(
                        id='title', 
                        children=app_config.params['title'], 
                        style=style_text_box
                    )]
            ),
            html.Div(
                className="row", 
                children=[
                    create_div_landscapes(point_names, feat_names, more_colorvars), 
                    create_div_sidepanels(point_names, feat_names, more_colorvars)
                ]
            ), 
            create_div_cosmetic_panel(), 
            html.Div([ html.Pre(id='test-select-data', style={ 'color': app_config.params['font_color'], 'overflowX': 'scroll' } ) ]),     # For testing purposes only!
            html.Div(
                className='row', 
                children=[ 
                    dcc.Markdown(
                        """Queries? Requests? Contact [Akshay Balsubramani](abalsubr@stanford.edu). """ 
                        + """Source [repository](https://github.com/kundajelab/coessentiality-browser)."""
                        )], 
                style={
                    'textAlign': 'center', 
                    'color': app_config.params['font_color'], 
                    'padding-bottom': '10px'
                }
            ), 
            dcc.Store(
                id='stored-pointsets', 
                data={ '_current_selected_data': {} }    # Maintained as the short-term state of a point subset.
            ), 
            dcc.Store(
                id='stored-landscape-selected', 
                data={ }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-most-recently-highlighted', 
                data={ '_last_panel_highlighted': 'landscape' }, 
                modified_timestamp='0'
            )
        ],
        style={ 
            'width': '100vw', 
            'max-width': 'none'
        }
    )
