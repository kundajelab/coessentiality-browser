# -*- coding: utf-8 -*-
"""
Runner for an interactive browser / analyzer.
Dash application (Dash [Python] <- Plotly <- React.js <- D3.js)
Author: Akshay Balsubramani
"""

import base64, io, os, time, json
import numpy as np, scipy as sp, pandas as pd, dash
from dash.dependencies import Input, Output, State
import app_config, app_lib, building_block_divs
# import matplotlib, matplotlib.pyplot as plt, matplotlib.colors as colors
"""
For more on jobs that take a while: set up workers https://github.com/WileyIntelligentSolutions/wiley-boilerplate-dash-app
"""


    

# =========================================================
# ================== Initialize Dash app ==================
# =========================================================


# Load gene embedded coordinates.
plot_data_df = pd.read_csv(app_config.params['plot_data_df_path'][0], sep="\t", index_col=False)
# graph_adj = sp.sparse.load_npz(app_config.params['adj_mat_path'])
data_ess = pd.read_csv(app_config.params['raw_ess_data_path'], index_col=0, header=0, sep='\t')
data_ess = data_ess[data_ess.columns[:-4]]   # Only first 481 cols put through GLS, so isolate these

point_names = np.array(plot_data_df['gene_names'])
feat_names = data_ess.columns
cancer_types = data_ess.columns.str.split('_').str[1:].str.join(' ').str.capitalize().str.replace('Haematopoietic and lymphoid tissue', 'Hematopoietic/lymphoid')
additional_colorvars = []#app_config.params['additional_colorvars']

raw_data = data_ess.values

app = dash.Dash(__name__)    #, external_stylesheets=external_stylesheets)
if not app_config._DEPLOY_LOCALLY:
    app.config.update({'routes_pathname_prefix':'/coessentiality/', 'requests_pathname_prefix':'/coessentiality/'})

server=app.server
app.layout = building_block_divs.create_div_mainapp(
    point_names, 
    feat_names, 
    more_colorvars=additional_colorvars
)



# =====================================================
# ===================== Callbacks =====================
# =====================================================


# =========================== The callback logic ===========================

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
            toret[p['text']] = {
                'pointIndex': p['pointIndex'], 
                'curveNumber': p['curveNumber']
            }
        return toret
    else:
        return {}

    
def make_selected(stored_dict):
    toret = { 'range': None }
    toret['points'] = [
        {
            'pointIndex': stored_dict[k]['pointIndex'], 
            'curveNumber': stored_dict[k]['curveNumber'], 
            'text': k, 
            'pointNumber': stored_dict[k]['pointIndex']
        } for k in stored_dict
    ]
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
    color_var=app_config.params['default_color_var'], # Could be an array of continuous colors!
    continuous_color=app_config.params['continuous_color'], 
    colorscale=app_config.params['colorscale'], 
    selectedpoint_ids=[], 
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
            'arrowcolor': 'white', 
            'showarrow': True, 
            'arrowhead': 2, 'arrowwidth': 2, 'arrowsize': 2, 
            'ax': 0, 
            'ay': -50 
        })
    toret = app_lib.build_main_scatter(
        data_df, 
        color_var, 
        colorscale, 
        annots=annots, 
        selected_point_ids=selectedpoint_ids, 
        bg_marker_size=bg_marker_size, 
        marker_size=marker_size
    )
    return toret


def run_update_main_heatmap(
    view_option, 
    subset_store, 
    landscape_scatter_fig, 
    point_names_to_use, 
    raw_data_to_use, 
    num_points_to_sample=10000
):
    pointIDs_to_display = list(subset_store['_current_selected_data'].keys())
    # Subsample down to a max #points, for smoothly interactive heatmap display.
    if len(pointIDs_to_display) > num_points_to_sample:
        pointIDs_to_display = np.random.choice(pointIDs_to_display, num_points_to_sample, replace=False)
    point_ndces_to_display = np.isin(point_names_to_use, pointIDs_to_display)
    subset_raw_data = raw_data_to_use[point_ndces_to_display, :]
    if sp.sparse.issparse(raw_data_to_use):
        subset_raw_data = subset_raw_data.toarray()
    subset_point_names = point_names_to_use[point_ndces_to_display]
    cocluster_mode = (view_option == 'Cocluster')
    hm_fig = app_lib.display_heatmap_cb(
        subset_raw_data, 
        feat_names, 
        subset_point_names, 
        landscape_scatter_fig, 
        cocluster_mode
    )
    return hm_fig


def run_update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points,      # Selected points annotated
    subset_store,       # Store of selected point subsets.
    highlight_selected_points, 
    data_df, 
    point_names, 
    raw_data_to_use, 
    bg_marker_size, 
    marker_size, 
    highlighted_points=[]
):
    pointIDs_to_select = highlighted_points if (len(highlighted_points) > 0) else list(subset_store['_current_selected_data'].keys())
    if annotated_points is None:
        annotated_points = []
    if (len(annotated_points) == 0) and highlight_selected_points:
        annotated_points = pointIDs_to_select
    absc_arr = data_df[app_config.params['display_coordinates']['x']]
    ordi_arr = data_df[app_config.params['display_coordinates']['y']]
    # Check if a continuous feature is chosen to be plotted.
    if ((color_scheme != app_config.params['default_color_var']) and 
        (color_scheme not in additional_colorvars) and 
        (len(color_scheme) > 0)
       ):
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
            continuous_color=True, 
            colorscale=app_config.params['colorscale_continuous'], 
            selectedpoint_ids=pointIDs_to_select, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            bg_marker_size=bg_marker_size, 
            marker_size=marker_size
        )
    else:    # color_scheme is a col ID indexing a discrete column.
        return highlight_landscape_func(
            annotated_points, 
            data_df, 
            point_names, 
            color_var=color_scheme, 
            continuous_color=False, 
            colorscale=app_config.params['colorscale_discrete'], 
            selectedpoint_ids=pointIDs_to_select, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            bg_marker_size=bg_marker_size, 
            marker_size=marker_size
        )

    
def parse_upload_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    return json.loads(decoded)
    if 'csv' in filename:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    elif 'xls' in filename:
        df = pd.read_excel(io.BytesIO(decoded))



# ==============================================================================
# =========================== Callbacks with headers ===========================
# ==============================================================================


@app.callback(
    Output('test-select-data', 'children'),
    [Input('toggle-debug-panels', 'values'), 
     Input('landscape-plot', 'selectedData'), 
     Input('stored-pointsets', 'data'), 
     Input('main-heatmap', 'selectedData'), 
     Input('stored-most-recently-highlighted', 'data')]
)
def display_test(
    debug_options, 
    sel_data, 
    data_store, 
    hmsel_data, 
    hlight_store
):
    see_hm = "0" if hmsel_data is None else str(len(hmsel_data['points']))
    see_sel = "0" if sel_data is None else str(len(sel_data['points']))
    see_hlight = "{}\t{}".format(hlight_store['_last_panel_highlighted'], len(hlight_store.keys()) - 1)
    toret = ""
    for setname in data_store:
        toret = toret + "{}\t{}\n".format(len(data_store[setname]), setname)
    if 'debug-panel' in debug_options:
        return "***STORED SELECTED DATA***\n{}\n***Landscape SELECTED DATA***\n{}\n***Heatmap SELECTED DATA***\n{}\n***Most recently used panel:\t{}".format(
            toret, 
            see_sel, 
            see_hm, 
            see_hlight
        )
    else:
        return ""


# Link currently selected data to output of gene set selector, so it can be picked too.
@app.callback(
    Output('goenrich-panel', 'figure'),
    [Input('geneset-select', 'value'), 
     Input('select-topk-goterms', 'n_submit'),
     Input('select-topk-goterms', 'n_blur'), 
     Input('stored-pointsets', 'data')],
    [State('select-topk-goterms', 'value')]
)
def display_goenrich_panel(selected_genes, dummy1, dummy2, subset_store, topk):
    if len(selected_genes) == 0:
        selected_genes = list(subset_store['_current_selected_data'].keys())
    return app_lib.display_goenrich_panel_func(selected_genes, topk=int(topk))


# https://community.plot.ly/t/download-raw-data/4700/7
@app.callback(
    Output('download-set-link', 'href'),
    [Input('stored-pointsets', 'data')]
)
def save_selection(subset_store):
    save_contents = json.dumps(subset_store['_current_selected_data'], indent=4)
    return "data:text/json;charset=utf-8," + save_contents


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
    Output('display-bg-marker-size-factor', 'children'), 
    [Input('slider-bg-marker-size-factor', 'value')]
)
def update_bgmarker_size(bg_size):
    return 'Unannotated marker size: {}'.format(bg_size)


@app.callback(
    Output('display-marker-size-factor', 'children'), 
    [Input('slider-marker-size-factor', 'value')]
)
def update_marker_size(marker_size):
    return 'Marker size: {}'.format(marker_size)


@app.callback(
    Output('stored-landscape-selected', 'data'), 
    [Input('landscape-plot', 'selectedData')]
)
def update_stored_landscape_data(landscape_data):
    return make_store_points(landscape_data)


@app.callback(
    Output('landscape-plot', 'selectedData'), 
    [Input('stored-pointsets', 'data')]
)
def update_selected_landscape_data(subset_store):
    return make_selected(subset_store['_current_selected_data'])


# The following determines what the most recently called auxiliary panel dictated to highlight.
@app.callback(
    Output('stored-most-recently-highlighted', 'data'), 
    [Input('stored-landscape-selected', 'modified_timestamp'), 
     Input('hm-highlight-button', 'n_clicks_timestamp')], 
    [State('landscape-plot', 'selectedData'), 
     State('main-heatmap', 'selectedData')]
)
def update_most_recently_highlighted(
    landscape_time, 
    hm_time, 
    landscape_selection, 
    hm_selection
):
    recent_times = np.array([
        int(landscape_time), int(hm_time), 0
    ])
    most_recent_time = max(recent_times)
    if most_recent_time == int(landscape_time):
        rec_panel = 'landscape'
        data_selection = landscape_selection
    elif most_recent_time == int(hm_time):
        rec_panel = 'heatmap'
        data_selection = hm_selection
    # TODO add more aux panels here. keys should be _last_panel_highlighted as well as cell IDs in highlighted subset.
    toret = make_store_points(data_selection)
    toret['_last_panel_highlighted'] = rec_panel
    return toret


import dash_core_components as dcc
import plotly.figure_factory as ff
@app.callback(
    Output('hm-feat-control', 'children'), 
    [Input('toggle-hm-feat-panels', 'values')]
)
def update_hm_control_panel(panel_list):
    graphs = []
    cell_cluster_list = np.array([])  # List of cluster IDs for resp. cells
    cell_color_list = np.array([])
    if 'bars' in panel_list:
        graphs.append(
            app_lib.generate_percluster_viz(raw_data, cell_cluster_list, cell_color_list))
    if 'dendrogram' in panel_list:
        X = np.random.rand(15, 15)
        dendro = ff.create_dendrogram(X)
        graphs.append(dcc.Graph(figure=dendro))
    return graphs


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


"""
Updates the stored dictionary of saved subsets. 
Contains control logic for subset selection and storage.
"""
@app.callback(
    Output('stored-pointsets', 'data'), 
    [Input('store-status', 'values'), 
     Input('list-pointsets', 'value'), 
     Input('stored-most-recently-highlighted', 'data'), 
     Input('upload-pointsets', 'contents')],
    [State('upload-pointsets', 'filename'), 
     State('load-status', 'values'), 
     State('pointset-name', 'value'), 
     State('stored-pointsets', 'data'), 
     State('toggle-hm-zoom', 'values')]
)
def update_subset_storage(
    store_status, 
    selected_subsetIDs, 
    aux_highlighted, 
    file_contents, 
    file_paths, 
    load_status, 
    newset_name, 
    subset_store, 
    zoom_hm_selected
):
    new_sets_dict = {} if subset_store is None else subset_store
    if 'load' in load_status:    # Update _current_selected_data by loading subsets.
        new_sets_dict['_current_selected_data'] = union_of_selections(selected_subsetIDs, subset_store)
    else:   # Update _current_selected_data from the main plot / heatmap.
        # Logic to display points as selected from an auxplot (most recently used for selection). 
        # A small subset selected_heatmap_points should not change selected_landscape_points, but should change _current_selected_data
        last_hlight_panel = aux_highlighted.pop('_last_panel_highlighted', None)
        print(len(aux_highlighted))
        if last_hlight_panel == 'landscape':
            new_sets_dict['_current_selected_data'] = aux_highlighted
        elif last_hlight_panel == 'heatmap':
            if 'on' in zoom_hm_selected:
                new_sets_dict['_current_selected_data'] = aux_highlighted
    # Store current selected data as a new set if in that mode.
    if 'store' in store_status:
        if ((newset_name is not None) and 
            (newset_name not in new_sets_dict) and 
            (newset_name != '')):
            new_sets_dict[newset_name] = new_sets_dict['_current_selected_data']
    # Load a bunch of cell sets with names equal to their filenames.
    if file_contents is not None and len(file_contents) > 0:
        for contents, fname in zip(file_contents, file_paths):   # fname here is a relative (NOT an absolute) file path
            fname_root = fname.split('/')[-1].split('.')[0]
            new_sets_dict[fname_root] = parse_upload_contents(contents, fname)
    return new_sets_dict


"""
Update the main heatmap.
"""
@app.callback(
    Output('main-heatmap', 'figure'),
    [Input('main-heatmap-roworder', 'value'), 
     Input('stored-pointsets', 'data')], 
    [State('landscape-plot', 'figure')]
)
def update_main_heatmap(
    view_option, 
    subset_store, 
    landscape_scatter_fig, 
    num_points_to_sample=10000
):
    return run_update_main_heatmap(
        view_option, 
        subset_store, 
        landscape_scatter_fig, 
        point_names, 
        raw_data, 
        num_points_to_sample=num_points_to_sample
    )


"""
Update the main graph panel.
"""
@app.callback(
    Output('landscape-plot', 'figure'), 
    [Input('landscape_color', 'value'), 
     Input('points_annot', 'value'), 
     Input('stored-pointsets', 'data'), 
     Input('sourcedata-select', 'value'), 
     Input('toggle-landscape-whiten', 'values'), 
     Input('stored-most-recently-highlighted', 'data'), 
     Input('slider-bg-marker-size-factor', 'value'), 
     Input('slider-marker-size-factor', 'value')]
)
def update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points,      # Selected points annotated
    subset_store,          # Store of selected point subsets.
    sourcedata_select, 
    whiten_hm_selected, 
    aux_highlighted, 
    bg_marker_size, 
    marker_size
):
    print('Color scheme: {}'.format(color_scheme))
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
    recently_highlighted = [x for x in aux_highlighted.keys()]
    recently_highlighted.remove('_last_panel_highlighted')
    return run_update_landscape(
        color_scheme, 
        annotated_points, 
        subset_store, 
        'on' in whiten_hm_selected, 
        data_df, 
        point_names, 
        raw_data, 
        bg_marker_size, 
        marker_size, 
        highlighted_points=recently_highlighted    # List of most recently highlighted cells.
    )



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    app.run_server(port=8053, debug=True)
