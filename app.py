# -*- coding: utf-8 -*-
"""
Runner for an interactive browser / analyzer.
Dash application (Dash [Python] <- Plotly <- React.js <- D3.js)
Author: Akshay Balsubramani
"""

import base64, io, os, time, json
import numpy as np, scipy as sp, pandas as pd, dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import app_config, app_lib, building_block_divs




# =========================================================
# ================== Initialize Dash app ==================
# =========================================================

# Load gene embedded coordinates.
plot_data_df = pd.read_csv(app_config.params['plot_data_df_path'][0], sep="\t", index_col=False)
plot_data_df['Colors'] = 'Unannotated'

# graph_adj = sp.sparse.load_npz(app_config.params['adj_mat_path'])
data_ess = pd.read_csv(app_config.params['raw_ess_data_path'], index_col=0, header=0, sep='\t')
data_ess = data_ess[data_ess.columns[:-4]]   # Only first 481 cols put through GLS, so isolate these

point_names = np.array(plot_data_df['gene_names'])
# full_gene_ensIDs = np.array(pd.read_csv(app_config.params['gene_ensID_path'], sep="\t", header=None, index_col=False)).flatten()
feat_names = data_ess.columns
cancer_types = data_ess.columns.str.split('_').str[1:].str.join(' ').str.capitalize().str.replace('Haematopoietic and lymphoid tissue', 'Hematopoietic/lymphoid')

mutation_data = pd.read_csv(app_config.params['mutation_arr_path'], sep='\t').values
# shRNA_df = pd.read_csv(app_config.params['shRNA_data_path'], sep=",", index_col=0)
expr_data = pd.read_csv(app_config.params['expression_arr_path'], sep='\t', header=0)

# go_termIDs = np.load(app_config.params['gotermIDs_path'])
# go_termnames = np.load(app_config.params['gotermnames_path'])

ctypes = [x for x in np.unique(cancer_types)]
colorlist = app_config.cmap_celltypes
cell_line_colordict = dict([x for x in zip(ctypes, colorlist[0:len(ctypes)])])

raw_data = data_ess.values


# ====================================================
# ===================== Main app =====================
# ====================================================

meta_tags=[
    # A description of the app, used by e.g.
    # search engines when displaying search results.
    {
        'name': 'description',
        'content': 'My description'
    },
    # A tag that tells Internet Explorer (IE)
    # to use the latest renderer version available
    # to that browser (e.g. Edge)
    {
        'http-equiv': 'X-UA-Compatible',
        'content': 'IE=edge'
    },
    # A tag that tells the browser not to scale
    # desktop widths to fit mobile screens.
    # Sets the width of the viewport (browser)
    # to the width of the device, and the zoom level
    # (initial scale) to 1.
    #
    # Necessary for "true" mobile support.
    {
      'name': 'viewport',
      'content': 'width=device-width, initial-scale=1.0'
    }
]

app = dash.Dash(__name__, meta_tags=meta_tags)    #, external_stylesheets=external_stylesheets)

if not app_config._DEPLOY_LOCALLY:
    app.config.update({'routes_pathname_prefix':'/coessentiality/', 'requests_pathname_prefix':'/coessentiality/'})

server=app.server
app.layout = building_block_divs.create_div_mainapp(
    point_names, 
    feat_names, 
    ctypes, 
    app.get_asset_url('upload.png'), 
    app.get_asset_url('download.png')
)
app.title = 'Co-essentiality browser'

# =====================================================
# ===================== Callbacks =====================
# =====================================================

# All the callback logic 

"""
Utility function to take union of a list of selected subsets that have been stored.
Returns: dictionary. Keys: point IDs in the union. Values: {'pointIndex': p, 'curveNumber': v }.
"""
def union_of_selections(selected_subsets, subset_store):
    points_so_far = np.array([])
    dict_to_return = {}
    if selected_subsets is not None:
        for v in selected_subsets:
            new_point_ids = np.setdiff1d(list(subset_store[v].keys()), points_so_far)
            points_so_far = np.union1d(points_so_far, new_point_ids)
            for c in new_point_ids:
                dict_to_return[c] = subset_store[v][c]
    return dict_to_return


def make_store_points(selectedData_points):
    if (selectedData_points is not None) and ('points' in selectedData_points):
        toret = {}
        for p in selectedData_points['points']:
            pt_txt = p['text'].split('<br>')[0]
            toret[pt_txt] = {}
        return toret
    else:
        return {}

    
def make_selected(stored_dict):
    toret = { 'range': None }
    toret['points'] = [ { 'text': k } for k in stored_dict ]
    return toret


"""
Update the main graph panel with selected points annotated, using the given dataset.
"""
def highlight_landscape_func(
    annotated_points, 
    data_df, 
    point_names_to_use, 
    style_selected=building_block_divs.style_selected, 
    color_var=app_config.params['default_color_var'], # Could be an array of continuous colors!
    colorscale=app_config.params['colorscale'], 
    selectedpoint_ids=[], 
    highlight_selected=False, 
    absc_arr=None, 
    ordi_arr=None
):
    annots = []
    looked_up_ndces = np.where(np.in1d(point_names_to_use, annotated_points))[0]
    for point_ndx in looked_up_ndces:
        absc = absc_arr[point_ndx]
        ordi = ordi_arr[point_ndx]
        cname = point_names_to_use[point_ndx]
        annots.append({
            'x': absc, 
            'y': ordi,
            'xref': 'x', 'yref': 'y', 
            # 'text': '<b>Cell {}</b>'.format(cname), 
            'font': { 
                'color': 'white', 
                'size': 15 
            }, 
            'arrowcolor': '#ff69b4', 
            'showarrow': True, 
            'arrowhead': 2, 'arrowwidth': 2, 'arrowsize': 2, 
            'ax': 0, 
            'ay': -50 
        })
    toret = app_lib.build_main_scatter(
        data_df, 
        color_var, 
        colorscale, 
        highlight=highlight_selected, 
        annots=annots, 
        selected_point_ids=selectedpoint_ids, 
        style_selected=style_selected, 
        point_names=point_names
    )
    return toret


def run_update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points,      # Selected points annotated
    subset_store,       # Store of selected point subsets.
    data_df, 
    point_names, 
    raw_data_to_use, 
    style_selected, 
    highlighted_points=[], 
    aggregate_tissue=False
):
    pointIDs_to_select = highlighted_points if (len(highlighted_points) > 0) else list(subset_store.keys())
    if annotated_points is None:
        annotated_points = []
    absc_arr = data_df[app_config.params['display_coordinates']['x']]
    ordi_arr = data_df[app_config.params['display_coordinates']['y']]
    # Check if a continuous feature is chosen to be plotted.
    if ((color_scheme is not None) and 
        (color_scheme != app_config.params['default_color_var']) and 
        (len(color_scheme) > 0)
       ):
        if aggregate_tissue:
            cell_line_ndces = np.where(cancer_types == color_scheme)[0]
            new_colors = np.array(raw_data_to_use[:, cell_line_ndces].mean(axis=1))
        else:
            if not isinstance(color_scheme, (list, tuple)):
                color_scheme = [color_scheme]
            feat_ndces = np.isin(feat_names, color_scheme)
            # If there are multiple continuous features given, plot their mean; otherwise the mean's a no-op anyway so it's safe to use.
            if sp.sparse.issparse(raw_data_to_use):
                new_colors = np.squeeze(np.array(raw_data_to_use[:, feat_ndces].mean(axis=1)))
            else:
                new_colors = np.mean(raw_data_to_use[:, feat_ndces], axis=1)
        return highlight_landscape_func(
            annotated_points, 
            data_df, 
            point_names, 
            color_var=new_colors, 
            colorscale=app_config.params['colorscale_continuous'], 
            selectedpoint_ids=pointIDs_to_select, 
            highlight_selected=True, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            style_selected=style_selected
        )
    else:    # color_scheme is a col ID indexing a discrete column.
        return highlight_landscape_func(
            annotated_points, 
            data_df, 
            point_names, 
            color_var=color_scheme, 
            colorscale=app_config.params['colorscale_discrete'], 
            selectedpoint_ids=pointIDs_to_select, 
            highlight_selected=True, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            style_selected=style_selected
        )

    
def parse_upload_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    return decoded.decode('utf-8').splitlines()
    if 'json' in filename:
        return json.loads(decoded)
    elif 'csv' in filename:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    elif 'xls' in filename:
        df = pd.read_excel(io.BytesIO(decoded))



# ================= Callbacks with headers ===================

# @app.callback(
#     Output('test-select-data', 'children'),
#     [Input('stored-panel-settings', 'data'), 
#      Input('stored-landscape-selected', 'data'), 
#      Input('stored-pointsets', 'data')]
# )
# def display_test(
#     panel_data, 
#     sel_data, 
#     data_store
# ):
#     see_sel = "0" if sel_data is None else str(len(sel_data))
#     toret = ""
#     for setname in data_store:
#         toret = toret + "{}\t{}\n".format(str(len(data_store[setname])), setname)
#     if panel_data['debug_panel']:
#         return "***STORED SELECTED DATA***\n{}\n***Landscape SELECTED DATA***\n{}".format(
#             toret, 
#             see_sel
#         )
#     else:
#         return ""


# https://community.plot.ly/t/download-raw-data/4700/7
@app.callback(
    Output('download-set-link', 'href'),
    [Input('landscape-plot', 'selectedData')]
)
def save_selection(landscape_data):
    subset_store = make_store_points(landscape_data)
    save_contents = '\n'.join(list(subset_store.keys()))
    return "data:text/csv;charset=utf-8," + save_contents


# Update dialogs.
@app.callback(
    Output('num-selected-counter', 'children'), 
    [Input('landscape-plot', 'selectedData')]
)
def update_numselected_counter(
    landscape_data
):
    subset_store = make_store_points(landscape_data)
    num_selected = len(subset_store)
    return '# selected: {}'.format(num_selected)


# Link currently selected data to output of gene set selector, so it can be picked too.
@app.callback(
    Output('goenrich-panel', 'figure'),
    [Input('hm-button', 'n_clicks'), 
     Input('select-topk-goterms', 'n_submit'),
     Input('select-topk-goterms', 'n_blur')],
    [State('landscape-plot', 'selectedData'), 
     State('select-topk-goterms', 'value')]
)
def display_goenrich_panel(
    button_numclicks, dummy1, dummy2, landscape_data, topk
):
    sel_hlight = make_store_points(landscape_data)
    selected_genes = [x for x in sel_hlight.keys()]
    return app_lib.display_goenrich_panel_func(selected_genes, topk=int(topk))


@app.callback(
    Output('select-geneview', 'options'), 
    [Input('select-hm-dataset', 'value')]
)
def update_geneview_options(geneview_dataset):
    gene_names = point_names
    if geneview_dataset == 'Mutation':
        gene_names = np.intersect1d(mutation_data[:,0], point_names)
    elif geneview_dataset == 'Expression':
        gene_names = np.intersect1d(expr_data.values[:,0], point_names)
    return [ {'value': gn, 'label': gn} for gn in gene_names ]


"""
Update the main heatmap.
"""
@app.callback(
    [Output('main-heatmap', 'figure'), 
     Output('main-heatmap', 'selectedData')],
    [Input('hm-button', 'n_clicks'), 
     Input('select-hm-dataset', 'value'), 
     Input('select-geneview', 'value')], 
    [State('landscape-plot', 'selectedData'), 
     State('cell-line-lookup', 'value'), 
     State('landscape-plot', 'figure')]
)
def update_main_heatmap(
    button_numclicks, 
    geneview_mode, 
    geneview_gene, 
    landscape_data, 
    landscape_color, 
    landscape_scatter_fig, 
    num_points_to_sample=10000
):
    geneview_celllines = None
    if geneview_mode == 'Mutation':
        geneview_data = mutation_data
    elif geneview_mode == 'Expression':
        geneview_data = expr_data.values
        geneview_celllines = expr_data.columns
    subset_store = make_store_points(landscape_data)
    col_names = feat_names
    point_names_to_use = point_names
    raw_data_to_use = raw_data
    
    pointIDs_to_display = list(subset_store.keys())
    # Subsample down to a max #points, for smoothly interactive heatmap display.
    if len(pointIDs_to_display) > num_points_to_sample:
        pointIDs_to_display = np.random.choice(pointIDs_to_display, num_points_to_sample, replace=False)
    # If any points (genes) are selected but not in the heatmap, they won't be displayed.
    point_ndces_to_display = np.isin(point_names_to_use, pointIDs_to_display)
    subset_point_names = point_names_to_use[point_ndces_to_display]
    subset_raw_data = raw_data_to_use[point_ndces_to_display, :]
    if sp.sparse.issparse(raw_data_to_use):
        subset_raw_data = subset_raw_data.toarray()
    cocluster_mode = False
    hm_fig = app_lib.display_heatmap_cb(
        subset_raw_data, 
        col_names, 
        subset_point_names, 
        landscape_scatter_fig, 
        cocluster_mode, 
        geneview_mode=geneview_mode, 
        geneview_gene=geneview_gene, 
        geneview_data=geneview_data, 
        geneview_celllines=geneview_celllines, 
        feat_group_names=cancer_types,
        feat_colordict=cell_line_colordict, 
        show_legend=False
    )
    return hm_fig, landscape_data#None


@app.callback(
    Output('landscape-plot', 'selectedData'), 
    [Input('main-heatmap', 'selectedData')]
)
def update_stored_heatmap_data(hm_selection):
    return hm_selection#make_store_points(hm_selection)


"""
Update the main graph panel.
"""
@app.callback(
    Output('landscape-plot', 'figure'), 
    [Input('cell-line-lookup', 'value'), 
     Input('tissue-type-lookup', 'value'), 
     Input('gene-selection', 'value'), 
     Input('landscape-plot', 'selectedData')]
)
def update_landscape(
    cell_line_color,          # Feature(s) selected to plot as color.
    selected_tissue_color,       # Selected tissue type (group of cell lines)
    annotated_points,      # Selected points annotated
    landscape_data
):
    subset_store = make_store_points(landscape_data)
    dataset_names = app_config.params['dataset_options']
    ndx_selected = 0 # dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
    data_df['Colors'] = 'Unannotated'
    style_selected = building_block_divs.style_selected
    if (selected_tissue_color is not None) and len(selected_tissue_color) > 0:
        color_scheme = selected_tissue_color
        aggregate_tissue = True
    else:
        color_scheme = cell_line_color if cell_line_color is not None else app_config.params['default_color_var']
        aggregate_tissue = False
    # print (selected_tissue_color, cell_line_color, color_scheme, aggregate_tissue)
    lscape = run_update_landscape(
        color_scheme, 
        annotated_points, 
        subset_store, 
        data_df, 
        point_names, 
        raw_data, 
        style_selected, 
        aggregate_tissue=aggregate_tissue
    )
    return lscape
        



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    #from waitress import serve
    #serve(server,host='0.0.0.0',port=5000)
    app.run_server(host='0.0.0.0', port=5000, debug=True)
