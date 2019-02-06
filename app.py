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

_DEPLOY_LOCALLY = True

# =========================================================
# ================== Initialize Dash app ==================
# =========================================================

if not _DEPLOY_LOCALLY:
    data_pfx = '/var/www/coessentiality-browser/'
else:
    data_pfx = '/Users/akshay/github/coessentiality-browser/'


# Load gene embedded coordinates.
plot_data_df = pd.read_csv(data_pfx + app_config.params['plot_data_df_path'], sep="\t", index_col=False)
# graph_adj = sp.sparse.load_npz(app_config.params['adj_mat_path'])
data_ess = pd.read_csv(data_pfx + 'data/essentiality.tsv.gz', index_col=0, header=0, sep='\t')
data_ess = data_ess[data_ess.columns[:-4]]   # Only first 481 cols put through GLS, so isolate these

point_names = np.array(plot_data_df['gene_names'])
feat_names = data_ess.columns
cancer_types = data_ess.columns.str.split('_').str[1:].str.join(' ').str.capitalize().str.replace('Haematopoietic and lymphoid tissue', 'Hematopoietic/lymphoid')
additional_colorvars = []#app_config.params['additional_colorvars']

raw_data = data_ess.values

app = dash.Dash(__name__)    #, external_stylesheets=external_stylesheets)
if not _DEPLOY_LOCALLY:
    app.config.update({'routes_pathname_prefix':'/coessentiality/', 'requests_pathname_prefix':'/coessentiality/'})

server=app.server
app.layout = building_block_divs.create_div_mainapp(
    point_names, 
    feat_names, 
    more_colorvars=additional_colorvars, 
    align_options_list=['Unaligned', 'Aligned']
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



"""
Update the main graph panel with selected points annotated, using the given dataset.
"""
def highlight_landscape_func(
    annotated_points, 
    color_var=app_config.params['default_color_var'], # Could be an array of continuous colors!
    continuous_color=app_config.params['continuous_color'], 
    colorscale=app_config.params['colorscale'], 
    point_names_to_use=None, 
    selectedpoint_ids=[], 
    absc_arr=None, 
    ordi_arr=None
):
    if point_names_to_use is None:
        point_names_to_use = point_names
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
        plot_data_df, 
        color_var, 
        colorscale, 
        annots=annots, 
        selected_point_ids=selectedpoint_ids
    )
    return toret


def run_update_main_heatmap(
    view_option, 
    subset_store, 
    landscape_scatter_fig, 
    num_points_to_sample=10000, 
    point_names_to_use=None, raw_data_to_use=None     # These two should match if either is provided
):
    if point_names_to_use is None:
        point_names_to_use = point_names
    if raw_data_to_use is None:
        raw_data_to_use = raw_data
    # Start by taking the currently selected points, and add points explicitly selected in landscape.
    pointIDs_to_display = list(subset_store['_current_selected_data'].keys())
#     if landscape_selected_data is not None:
#         these_selections = [ p['text'] for p in landscape_selected_data['points'] ]
#         pointIDs_to_display = np.union1d(pointIDs_to_display, these_selections)
    
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
    raw_data_to_use=None
):
    if raw_data_to_use is None:
        raw_data_to_use = raw_data
    pointIDs_to_select = list(subset_store['_current_selected_data'].keys())
    if annotated_points is None:
        annotated_points = []
    if (len(annotated_points) == 0) and 'highlight' in highlight_selected_points:
        annotated_points = pointIDs_to_select
    absc_arr = plot_data_df[app_config.params['display_coordinates']['x']]
    ordi_arr = plot_data_df[app_config.params['display_coordinates']['y']]
    
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
            color_var=new_colors, 
            continuous_color=True, 
            colorscale=app_config.params['colorscale_continuous'], 
            selectedpoint_ids=pointIDs_to_select, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr
        )
    else:    # color_scheme is a col ID indexing a discrete column.
        colorscale = app_config.params['colorscale_discrete']
        return highlight_landscape_func(
            annotated_points, 
            color_var=color_scheme, 
            continuous_color=False, 
            colorscale=colorscale, 
            selectedpoint_ids=pointIDs_to_select, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr
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


"""
# FOR TESTING ONLY
@app.callback(
    Output('test-select-data', 'children'),
    [Input('landscape-plot', 'selectedData'), 
     Input('stored-pointsets', 'data'), 
     Input('main-heatmap', 'selectedData')]
)
def display_test(
    sel_data, 
    data_store, 
    hmsel_data
):
    toret = ""
    see_hm = "0" if hmsel_data is None else str(len(hmsel_data['points']))
    see_sel = "0" if sel_data is None else str(len(sel_data['points']))
    for setname in data_store:
        toret = toret + "{}\t{}\n".format(len(data_store[setname]), setname)
    return "***STORED SELECTED DATA***\n{}\n***Landscape SELECTED DATA***\n{}\n***Heatmap SELECTED DATA***\n{}".format(
        toret, 
        see_sel, 
        see_hm
    )
"""

#,[State('select-topk-goterms', 'value')]
@app.callback(
    Output('goenrich-panel', 'figure'),
    [Input('stored-pointsets', 'data'),
     Input('select-topk-goterms', 'n_submit'),
     Input('select-topk-goterms', 'n_blur')],
    [State('select-topk-goterms', 'value')]
)
def display_goenrich_panel(subset_store, dummy1, dummy2, topk):
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
    toret = [] if stored_setlist is None else [ {'label': s, 'value': s} for s in stored_setlist.keys()]
    return toret


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
    [Input('store-button', 'n_clicks'), 
     Input('list-pointsets', 'value'), 
     Input('landscape-plot', 'selectedData'), 
     Input('upload-pointsets', 'contents')],
    [State('upload-pointsets', 'filename'), 
     State('load-status', 'values'), 
     State('pointset-name', 'value'), 
     State('stored-pointsets', 'data')]
)
def update_subset_storage(
    store_status, 
    selected_subsetIDs, 
    selected_landscape_points, 
    file_contents, 
    file_paths, 
    load_status, 
    newset_name, 
    subset_store
):
    new_sets_dict = {} if subset_store is None else subset_store
    if 'load' in load_status:    # Update _current_selected_data with union of loaded subsets if applicable.
        new_sets_dict['_current_selected_data'] = union_of_selections(selected_subsetIDs, subset_store)
    else:   # Update _current_selected_data with new selected data, from the main plot / heatmap.
        new_sets_dict['_current_selected_data'] = make_store_points(selected_landscape_points)
    # Finally store current selected data as a new set if in that mode.
    # if 'store' in store_status:
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


@app.callback(
    Output('landscape-plot', 'selectedData'), 
    [Input('main-heatmap', 'selectedData'), 
     Input('toggle-heatmap-selection', 'values')], 
    [State('landscape-plot', 'selectedData')]
)
def update_landscape_seldata(hm_selected, hm_override_status, old_ls_data):
    if 'hm_override' in hm_override_status:
        return hm_selected
    else:
        return old_ls_data


"""
Update the main heatmap.
"""
@app.callback(
    Output('main-heatmap', 'figure'),
    [Input('main-heatmap-roworder', 'value'), 
     Input('stored-pointsets', 'data'), 
     Input('landscape-plot', 'figure')]
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
     Input('toggle-heatmap-selection', 'values')]
)
def update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points,      # Selected points annotated
    subset_store,       # Store of selected point subsets.
    highlight_selected_points
):
    print('Color scheme: {}'.format(color_scheme))
    return run_update_landscape(
        color_scheme, annotated_points, subset_store, highlight_selected_points
    )



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    app.run_server(port=8053, debug=True)
