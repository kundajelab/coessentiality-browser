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
full_gene_ensIDs = np.array(pd.read_csv(app_config.params['gene_ensID_path'], sep="\t", header=None, index_col=False)).flatten()
feat_names = data_ess.columns
cancer_types = data_ess.columns.str.split('_').str[1:].str.join(' ').str.capitalize().str.replace('Haematopoietic and lymphoid tissue', 'Hematopoietic/lymphoid')

mutation_data = np.load(app_config.params['mutation_arr_path'], allow_pickle=True)
# shRNA_df = pd.read_csv(app_config.params['shRNA_data_path'], sep=",", index_col=0)

expr_data = np.load(app_config.params['expression_arr_path'], allow_pickle=True)
expr_cell_lines = np.load(app_config.params['expression_cell_lines_path'])

go_termIDs = np.load(app_config.params['gotermIDs_path'])
go_termnames = np.load(app_config.params['gotermnames_path'])

ctypes = [x for x in np.unique(cancer_types)]
colorlist = app_config.cmap_celltypes
cell_line_colordict = dict([x for x in zip(ctypes, colorlist[0:len(ctypes)])])

additional_colorvars = []# ['roarke_clusters']

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
    app.get_asset_url('download.png'), 
    full_gene_ensIDs, 
    more_colorvars=additional_colorvars
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
    bg_marker_size=app_config.params['bg_marker_size_factor'], 
    marker_size=app_config.params['marker_size_factor'], 
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
    if color_var is None:
        color_var = app_config.params['default_color_var']
    toret = app_lib.build_main_scatter(
        data_df, 
        color_var, 
        colorscale, 
        highlight=highlight_selected, 
        annots=annots, 
        selected_point_ids=selectedpoint_ids, 
        bg_marker_size=bg_marker_size, 
        marker_size=marker_size, 
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
    bg_marker_size, 
    marker_size, 
    style_selected, 
    highlighted_points=[], 
    aggregate_tissue=False
):
    pointIDs_to_select = highlighted_points if (len(highlighted_points) > 0) else list(subset_store['_current_selected_data'].keys())
    if annotated_points is None:
        annotated_points = []
    absc_arr = data_df[app_config.params['display_coordinates']['x']]
    ordi_arr = data_df[app_config.params['display_coordinates']['y']]
    # Check if a continuous feature is chosen to be plotted.
    if ((color_scheme is not None) and 
        (color_scheme != app_config.params['default_color_var']) and 
        (color_scheme not in additional_colorvars) and 
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
            bg_marker_size=bg_marker_size, 
            marker_size=marker_size, 
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
            bg_marker_size=bg_marker_size, 
            marker_size=marker_size, 
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

# Render selectable point subsets.
@app.callback(
    Output('list-pointsets', 'options'), 
    [Input('stored-pointsets', 'data')]
)
def update_subset_options(stored_setlist):
    toret = []
    blacklist = ['_current_subselected_data', '_current_selected_data']
    if stored_setlist is not None:
        for s in stored_setlist.keys():
            if s not in blacklist:
                toret.append({'label': s, 'value': s})
    return toret


@app.callback(
    Output('test-select-data', 'children'),
    [Input('stored-panel-settings', 'data'), 
     Input('stored-landscape-selected', 'data'), 
     Input('stored-pointsets', 'data'), 
     Input('stored-heatmap-selected', 'data')]
)
def display_test(
    panel_data, 
    sel_data, 
    data_store, 
    hm_store
):
    see_sel = "0" if sel_data is None else str(len(sel_data))
    see_hm = "0" if hm_store is None else str(len(hm_store))
    toret = ""
    for setname in data_store:
        toret = toret + "{}\t{}\n".format(str(len(data_store[setname])), setname)
    if panel_data['debug_panel']:
        return "***STORED SELECTED DATA***\n{}\n***Landscape SELECTED DATA***\n{}\n***Heatmap SELECTED DATA***\n{}".format(
            toret, 
            see_sel, 
            see_hm
        )
    else:
        return ""


# https://community.plot.ly/t/download-raw-data/4700/7
@app.callback(
    Output('download-set-link', 'href'),
    [Input('stored-pointsets', 'data')]
)
def save_selection(subset_store):
    save_contents = '\n'.join(list(subset_store['_current_selected_data'].keys()))
    return "data:text/csv;charset=utf-8," + save_contents
#     save_contents = json.dumps(subset_store['_current_selected_data'], indent=4)
#     return "data:text/json;charset=utf-8," + save_contents


import urllib
# Downloads the currently selected dataframe before extracting and exporting the desired columns.
@app.callback(
    Output('download-layout-link', 'href'), 
    [Input('sourcedata-select', 'value')]
)
def update_download_layout_link(sourcedata_select):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
    coords_to_load = list(app_config.params['display_coordinates'].values()) + ['gene_names']
    csvString = data_df[coords_to_load].to_csv(sep="\t", index=False,encoding='utf-8')
    return "data:text/csv;charset=utf-8," + csvString


# @app.callback(
#     Output('display-bg-marker-size-factor', 'children'), 
#     [Input('slider-bg-marker-size-factor', 'value')]
# )
# def update_bgmarker_size(bg_size):
#     return 'Unannotated marker size: {}'.format(bg_size)


# @app.callback(
#     Output('display-marker-size-factor', 'children'), 
#     [Input('slider-marker-size-factor', 'value')]
# )
# def update_marker_size(marker_size):
#     return 'Marker size: {}'.format(marker_size)


@app.callback(
    Output('stored-landscape-selected', 'data'), 
    [Input('landscape-plot', 'selectedData')]
)
def update_stored_landscape_data(landscape_data):
    return make_store_points(landscape_data)


@app.callback(
    Output('stored-heatmap-selected', 'data'), 
    [Input('main-heatmap', 'selectedData')]
)
def update_stored_heatmap_data(hm_selection):
    return make_store_points(hm_selection)


@app.callback(
    Output('landscape-plot', 'selectedData'), 
    [Input('stored-pointsets', 'data')]
)
def update_selected_landscape_data(subset_store):
    return make_selected(subset_store['_current_selected_data'])


# Handle lookups of GO terms and return a gene set.
@app.callback(
    Output('stored-goterm-lookup-results', 'data'), 
    [Input('goterm-lookup', 'value')]
)
def update_goterm_lookup(
    goterms_req	
):
    tmpl = [app_lib.get_genes_from_goterm(termID) for termID in goterms_req]
    if len(tmpl) == 0:
        return ""
    sel_genes = np.concatenate(tmpl)
    return list(np.unique(sel_genes))


# Updates the stored dictionary of boolean panel config variables.
@app.callback(
    [Output('stored-recently-highlighted', 'data'), 
     Output('goterm-lookup', 'value'), 
     Output('time-recently-highlighted', 'data')], 
    [Input('stored-landscape-selected', 'modified_timestamp'), 
     Input('stored-heatmap-selected', 'modified_timestamp'), 
     Input('stored-goterm-lookup-results', 'modified_timestamp')], 
    [State('stored-landscape-selected', 'data'), 
     State('stored-heatmap-selected', 'data'), 
     State('stored-goterm-lookup-results', 'data'), 
     State('goterm-lookup', 'value')]
)
def update_hlight_data_store(
    time_sel_landscape, 
    time_sel_heatmap, 
    time_sel_goterm, 
    sel_landscape, 
    sel_heatmap, 
    sel_goterm, 
    goterm_lookup
):
    times_list = [int(time_sel_landscape), int(time_sel_heatmap), int(time_sel_goterm)]
    stores_list = [sel_landscape, sel_heatmap, { x: {} for x in sel_goterm }]
    times_list.append(0)
    most_recent_time_index = np.argmax(times_list)
    if most_recent_time_index < len(times_list):
        dt = []
        if most_recent_time_index == 2:
            dt = goterm_lookup
        # print(stores_list[most_recent_time_index], dt, times_list[most_recent_time_index])
        return stores_list[most_recent_time_index], dt, times_list[most_recent_time_index]
    else:
        return {}, [], 0


@app.callback(
    Output('select-geneview', 'options'), 
    [Input('select-hm-dataset', 'value')]
)
def update_geneview_options(geneview_dataset):
    gene_names = point_names
    if geneview_dataset == 'Mutation':
        gene_names = np.intersect1d(mutation_data[:,0], point_names)
    elif geneview_dataset == 'Expression':
        gene_names = np.intersect1d(expr_data[:,0], point_names)
    return [ {'value': gn, 'label': gn} for gn in gene_names ]

"""
@app.callback(
    [Output('goterm-lookup', 'options'), 
     Output('load-go-db', 'style'), 
     Output('goterm-lookup', 'disabled')], 
    [Input('load-go-button', 'n_clicks')], 
    [State('goterm-lookup', 'options')]
)
def update_go_db(
    button_clicks, 
    options
):
    if (button_clicks == 0) or ((options is not None) and len(options) > 0):
        raise PreventUpdate
    else:
        return (
            [{'value': '{}'.format(go_termIDs[i]), 'label': '{}: \t{}'.format(go_termIDs[i], go_termnames[i])} for i in range(len(go_termIDs)) ], 
            {'display': 'none'}, 
            False
        )#'Search GO')
"""
# def update_goterm_options(search_val, val, options):
#     return [{'value': '{}'.format(go_termIDs[i]), 'label': '{}: \t{}'.format(go_termIDs[i], go_termnames[i])} for i in range(len(go_termIDs)) ]
    # Make sure that the set values are in the option list, else they will disappear from the shown select list, but still part of the `value`.
    # return [ o for o in options if search_val in o["label"] or o["value"] in (val or []) ]


# @app.callback(
#     Output('hm-future-panels', 'children'), 
#     [Input('toggle-future-panels', 'value')]
# )
# def update_future_panels(panel_list):
#     graphs = []
#     cell_cluster_list = np.array([])  # List of cluster IDs for resp. cells
#     cell_color_list = np.array([])
#     if 'diff_features' in panel_list:
#         graphs.append(building_block_divs.create_div_diff_features())
#     return graphs


# Update dialogs.
@app.callback(
    Output('num-selected-counter', 'children'), 
    [Input('stored-pointsets', 'data')]
)
def update_numselected_counter(
    subset_store
):
    num_selected = len(subset_store['_current_selected_data'])
    return '# selected: {}'.format(num_selected)


# Link currently selected data to output of gene set selector, so it can be picked too.
@app.callback(
    Output('goenrich-panel', 'figure'),
    [Input('stored-recently-highlighted', 'data'), 
     Input('select-topk-goterms', 'n_submit'),
     Input('select-topk-goterms', 'n_blur')],
    [State('select-topk-goterms', 'value')]
)
def display_goenrich_panel(
    sel_hlight, dummy1, dummy2, topk
):
    selected_genes = [x for x in sel_hlight.keys()]
    return app_lib.display_goenrich_panel_func(selected_genes, topk=int(topk))


"""
Updates the stored dictionary of saved subsets. 
Contains control logic for subset selection and storage.
"""
@app.callback(
    Output('stored-pointsets', 'data'), 
    [Input('list-pointsets', 'value'), 
     Input('upload-pointsets', 'contents'), 
     Input('stored-recently-highlighted', 'data')],
    [State('upload-button', 'n_clicks_timestamp'), 
     State('upload-pointsets', 'filename'), 
     State('time-recently-highlighted', 'data'), 
     State('stored-pointsets', 'data')]
)
def update_subset_storage(
    selected_subsetIDs, 
    file_contents, 
    selected_hlight_data, 
    upload_time, 
    file_paths, 
    time_hlight_data, 
    subset_store
):
    new_sets_dict = {} if subset_store is None else subset_store
    if selected_subsetIDs is not None and len(selected_subsetIDs) > 0:    # Update _current_selected_data by loading subsets.
        new_sets_dict['_current_selected_data'] = union_of_selections(selected_subsetIDs, subset_store)
    else:   # Update _current_selected_data from the main plot / heatmap.
        # Logic to display points as selected from an auxplot (most recently used for selection). 
        # A small subset selected_heatmap_points should not change selected_landscape_points, but should change _current_selected_data
        new_sets_dict['_current_selected_data'] = selected_hlight_data
    # Store current selected data as a new set if in that mode.
    # if ((newset_name is not None) and (newset_name not in new_sets_dict) and (newset_name != '')):
    #     new_sets_dict[newset_name] = new_sets_dict['_current_selected_data']
    # Load a bunch of cell sets with names equal to their filenames.
    if (time_hlight_data is not None) and (upload_time > time_hlight_data) and (file_contents is not None and len(file_contents) > 0):
        for contents, fname in zip(file_contents, file_paths):   # fname here is a relative (NOT an absolute) file path
            fname_root = fname.split('/')[-1].split('.')[0]
            # Now make and save a new subset.
            new_sets_dict[fname_root] = { x: {} for x in parse_upload_contents(contents, fname) }
        new_sets_dict['_current_selected_data'] = new_sets_dict[fname_root]
    return new_sets_dict


"""
Update the main heatmap.
"""
@app.callback(
    Output('main-heatmap', 'figure'),
    [Input('stored-landscape-selected', 'data'), 
     Input('toggle-hm-cols', 'value'), 
     Input('select-hm-dataset', 'value'), 
     Input('select-geneview', 'value'), 
     Input('cell-line-lookup', 'value')], 
    [State('landscape-plot', 'figure')]
)
def update_main_heatmap(
    subset_store, 
    hm_col_panel, 
    geneview_mode, 
    geneview_gene, 
    landscape_color, 
    landscape_scatter_fig, 
    num_points_to_sample=10000
):
    if geneview_mode == 'Mutation':
        geneview_data = mutation_data
    elif geneview_mode == 'Expression':
        geneview_data = expr_data
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
        feat_group_names=cancer_types,
        feat_colordict=cell_line_colordict, 
        show_legend=False
    )
    return hm_fig


"""
Ensure at most one cell line / tissue lookup selected.
"""
# @app.callback(
#     Output('cell-line-lookup', 'value'), 
#     [Input('tissue-type-lookup', 'value')], 
#     [State('cell-line-lookup', 'value')]
# )
# def update_cell_line(selected_tissue_color, cell_line_color):
#     if (selected_tissue_color is not None) and len(selected_tissue_color) > 0:
#         return None
#     else:
#         return cell_line_color

@app.callback(
    [Output('tissue-type-lookup', 'value'), 
     Output('cell-line-lookup', 'value')], 
    [Input('tissue-or-cell-line', 'data')], 
    [State('tissue-type-lookup', 'value'), 
     State('cell-line-lookup', 'value')]
)
def update_color_dropdowns(
    selected_color, 
    selected_tissue_color, 
    selected_cell_line_color
):
    print(selected_color)
    if selected_color == 'cell-line':
        return '', selected_cell_line_color
    else:
        return selected_tissue_color, app_config.params['default_color_var']


@app.callback(
    Output('tissue-or-cell-line', 'data'), 
    [Input('div-tissue-lookup', 'n_clicks_timestamp'), 
     Input('div-cell-line-lookup', 'n_clicks_timestamp')], 
    [State('tissue-or-cell-line', 'data')]
)
def update_cell_line_tissue(
    time_tissue, 
    time_cellline, 
    old_state
):
    if (time_tissue is not None) and (time_cellline is not None):
        if (time_tissue > time_cellline) and (old_state != 'tissue'):
            return 'tissue'
        elif (time_tissue < time_cellline) and (old_state != 'cell-line'):
            return 'cell-line'
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate


# @app.callback(
#     Output('tissue-type-lookup', 'value'), 
#     [Input('cell-line-lookup', 'value')], 
#     [State('tissue-type-lookup', 'value')]
# )
# def update_cell_line(cell_line_color, selected_tissue_color):
#     if (cell_line_color is not None) and len(cell_line_color) > 0:
#         return None
#     else:
#         return selected_tissue_color


"""
Update the main graph panel.
"""
@app.callback(
    Output('landscape-plot', 'figure'), 
    [Input('cell-line-lookup', 'value'), 
     Input('tissue-type-lookup', 'value'), 
     Input('points_annot', 'value'), 
     Input('stored-pointsets', 'data'), 
     Input('sourcedata-select', 'value'), 
     # Input('select-ppi', 'value'), 
     Input('slider-bg-marker-size-factor', 'value'), 
     Input('slider-marker-size-factor', 'value')]
)
def update_landscape(
    cell_line_color,          # Feature(s) selected to plot as color.
    selected_tissue_color,       # Selected tissue type (group of cell lines)
    annotated_points,      # Selected points annotated
    subset_store,          # Store of selected point subsets.
    sourcedata_select, 
    # ppi_selected, 
    bg_marker_size, 
    marker_size
):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
    data_df['Colors'] = 'Unannotated'
    style_selected = building_block_divs.style_selected
    if (selected_tissue_color is not None) and len(selected_tissue_color) > 0:
        color_scheme = selected_tissue_color
        aggregate_tissue = True
    else:
        color_scheme = cell_line_color
        aggregate_tissue = False
    lscape = run_update_landscape(
        color_scheme, 
        annotated_points, 
        subset_store, 
        data_df, 
        point_names, 
        raw_data, 
        bg_marker_size, 
        marker_size, 
        style_selected, 
        aggregate_tissue=aggregate_tissue
    )
    # Add landscape edges as necessary.  # NO, too slow!
    itime = time.time()
    ppi_selected = 'None'
    if ppi_selected != 'None':
        if ppi_selected == 'STRING (v11)':
            adj_mat = sp.sparse.load_npz(app_config.params['string_matrix_ascoess_path']).tocoo()
        if ppi_selected == 'Cheng et al.':
            adj_mat = sp.sparse.load_npz(app_config.params['cheng_matrix_path']).tocoo()
        # Map point id to (x, y) pair
        pointIDs_to_coords = {}
        for trace in lscape['data']:
            for i in range(len(trace['x'])):
                pointIDs_to_coords[trace['text'][i]] = (trace['x'][i], trace['y'][i])
        
        edges_x = []
        edges_y = []
        row_names = point_names[adj_mat.row]
        col_names = point_names[adj_mat.col]
        for j in range(len(row_names)):
            coords1 = pointIDs_to_coords[row_names[j]]
            coords2 = pointIDs_to_coords[col_names[j]]
            edges_x.append((coords1[0], coords2[0], None))
            edges_y.append((coords1[1], coords2[1], None))
        lscape['data'].append({ 
            'name': 'Data', 
            'x': edges_x, 
            'y': edges_y, 
            'line': { 'width': 0.5, 'color': 'white'},
            'hoverinfo': 'none', 
            'mode': 'lines', 
            # 'render_mode': 'svg', 
            'type': 'scattergl'
        })
    return lscape
        



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    from waitress import serve
    serve(app,host='localhost',port=8080)
    #app.run_server(port=8052, debug=True)
