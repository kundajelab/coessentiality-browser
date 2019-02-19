# Application-specific routines for working with (smallish matrices of) data.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd, diffmap as dm
import app_config, building_block_divs
import umap
from sklearn.decomposition import TruncatedSVD



# =======================================================
# ================== Utility functions ==================
# =======================================================

# Loads anndata from 3 files: {data, rowdata, coldata}
def read_as_anndata(path_pfx, transp=True):
    m = scipy.io.mmread(path_pfx + '_data.mtx')
    if transp:
        m = m.transpose()
    c_annot = pd.read_csv(path_pfx + '_mdata_points.csv', index_col=None, header=0)
    f_annot = pd.read_csv(path_pfx + '_mdata_feats.csv', index_col=None, header=0)
    return anndata.AnnData(X=m, obs=c_annot, var=f_annot)


# Writes anndata in the form of 3 dataframes: {data, rowdata, coldata}
def write_anndata(path_pfx, adata, transp=True):
    if transp:
        adata = adata.transpose()
    adata.obs.to_csv(path_pfx + '_mdata_points.csv', index=False)
    adata.var.to_csv(path_pfx + '_mdata_feats.csv', index=False)
    scipy.io.mmwrite(path_pfx + '_data.mtx', adata.X)


def rna_log1pnorm(datamat, scale_mrna_per_point=None):
    mrna_per_point = np.sum(datamat, axis=1)
    if scale_mrna_per_point is None:
        scale_mrna_per_point = 10000
        # np.percentile(mrna_per_point, q=90)
    scale_factors = scale_mrna_per_point/mrna_per_point
    points_normed = np.multiply(datamat.T, scale_factors).T
    return np.log(1 + points_normed)


# Analysis mini-workflow for RNA-seq data
def fastfilt_rna(fdata, min_points_per_gene=25, max_expr_qtile_per_gene=99.5, scale_mrna_per_point=10000.0):
    mrna_per_point = np.sum(fdata.values, axis=1)
    mrna_per_gene = np.sum(fdata.values, axis=0)
    
    min_mrna_per_point = np.percentile(mrna_per_point, q=2)
    max_mrna_per_gene = np.percentile(mrna_per_gene, q=max_expr_qtile_per_gene)
    npoints, ngenes = fdata.values.shape
    point_mask = np.ones(npoints, dtype=bool)
    gene_mask = np.ones(ngenes, dtype=bool)
    point_mask = point_mask & (np.sum(fdata.values, axis=1) >= min_mrna_per_point)
    gene_mask = gene_mask & (np.sum(fdata.values > 0, axis=0) >= min_points_per_gene) 
    gene_mask = gene_mask & (np.sum(fdata.values, axis=0) <= max_mrna_per_gene)
    print("Fraction of points kept:\t" + str(np.mean(point_mask)))
    print("Fraction of genes kept:\t" + str(np.mean(gene_mask)))
    new_data = fdata[point_mask][fdata.columns[gene_mask]]
    new_data.iloc[:,:] = rna_log1pnorm(new_data.values)
    return new_data


def quantile_norm(dtd):
    qtiles = np.zeros(len(dtd))
    nnz_ndces = np.nonzero(dtd)[0]
    qtiles[nnz_ndces] = sp.stats.rankdata(dtd[nnz_ndces])/len(nnz_ndces)
    return qtiles


def interesting_feat_ndces(fit_data, num_feats_todisplay=500):
    num_feats_todisplay = min(fit_data.shape[1], num_feats_todisplay)
    if ((fit_data is None) or 
        (np.prod(fit_data.shape) == 0)
       ):
        return np.arange(num_feats_todisplay)
    feat_ndces = np.argsort(np.std(fit_data, axis=0))[::-1][:num_feats_todisplay]
    return feat_ndces


def reviz_embed_data(fit_data, alg='UMAP'):
#     feat_ndces = interesting_feat_ndces(fit_data)
#     fit_data = fit_data[:, feat_ndces]
    if alg == 'UMAP':
        reducer = umap.UMAP(n_neighbors=15, n_components=2, random_state=42)
        embedding = reducer.fit_transform(fit_data)
    return embedding


def pls_align(graph_adj, batch_IDs, n_components=2, umap=False):
    if n_components > 2:
        umap = True
    toret = np.zeros((len(batch_IDs), n_components))
    batch1_ndces = np.where(batch_IDs == np.unique(batch_IDs)[0])[0]
    batch2_ndces = np.where(batch_IDs == np.unique(batch_IDs)[1])[0]
    gg1 = graph_adj[batch1_ndces][:, batch1_ndces]
    gg2 = graph_adj[batch2_ndces][:, batch2_ndces]
    dmp_all_heat_1 = dm.diffmap_proj(gg1, n_comps=50, n_dims=20, min_energy_frac=0.9, embed_type='diffmap', return_eigvals=False)
    dmp_all_heat_2 = dm.diffmap_proj(gg2, n_comps=50, n_dims=20, min_energy_frac=0.9, embed_type='diffmap', return_eigvals=False)
    # dmp_all_naive, eigvals = dm.diffmap_proj(gg1, t=0, n_comps=reduced_dims, embed_type='naive', return_eigvals=True)
    cm = np.dot(dmp_all_heat_1, dmp_all_heat_2.T)
    tsvd = TruncatedSVD(n_components=n_components, random_state=42)
    toret[batch1_ndces, :] = tsvd.fit_transform(cm)
    toret[batch2_ndces, :] = tsvd.fit_transform(cm.T)
    if umap:
        pass
    return toret



# =========================================================
# =================== Main scatter plot ===================
# =========================================================

"""
(Data, layout) for the main graph panel.
Color_var is either a field of the plotting df, or a numpy array.
"""


# TODO: Add legend groups as applicable, to bunch colors within a group

# Here selected_point_ids is a list of unique string IDs of points. 
def traces_scatter(
    data_df, 
    color_var, 
    colorscale, 
    selected_point_ids, 
    bg_marker_size=app_config.params['bg_marker_size_factor'], 
    marker_size=app_config.params['marker_size_factor'], 
    style_selected=building_block_divs.style_selected
):
    traces_list = []
    display_ndces = app_config.params['display_coordinates']
    cumu_color_dict = {}
    # Check to see if color_var is continuous or discrete and plot points accordingly
    if isinstance(color_var, (list, tuple, np.ndarray)):     # Color_var is an array, not a col index.
        continuous_color_var = color_var
        point_names = list(data_df['gene_names'])
        spoints = np.where(np.isin(point_names, selected_point_ids))[0]
        colorbar_title = app_config.params['hm_colorvar_name']
        if app_config.params['qnorm_plot']:
            continuous_color_var = quantile_norm(continuous_color_var)
            colorbar_title = 'Percentile'
        pt_text = ["{}<br>Quantile: {}".format(point_names[i], round(continuous_color_var[i], 3)) for i in range(len(point_names))]
        traces_list.append({ 
            'name': 'Data', 
            'x': data_df[display_ndces['x']], 
            'y': data_df[display_ndces['y']], 
            'selectedpoints': spoints, 
            'hoverinfo': 'text', 
            'text': point_names, 
            'mode': 'markers', 
            'marker': {
                'size': bg_marker_size, 
                'opacity': app_config.params['bg_marker_opacity_factor'], 
                'symbol': 'circle', 
                'showscale': True, 
                'colorbar': {
                    'len': 0.3, 
                    'thickness': 20, 
                    'xanchor': 'right', 
                    'yanchor': 'top', 
                    'title': colorbar_title,
                    'titleside': 'top',
                    'ticks': 'outside', 
                    'titlefont': building_block_divs.colorbar_font_macro, 
                    'tickfont': building_block_divs.colorbar_font_macro
                }, 
                'color': continuous_color_var, 
                'colorscale': colorscale
            }, 
            'selected': style_selected, 
            'type': 'scattergl'
        })
    else:    # Categorical color scheme, one trace per color
        cnt = 0
        print('Colorscale length: {}'.format(len(colorscale)))
        for idx, val in data_df.groupby(color_var):
            point_ids_this_trace = list(val['gene_names'])
            spoint_ndces_this_trace = np.where(np.isin(point_ids_this_trace, selected_point_ids))[0]
            if idx not in cumu_color_dict:
                trace_color = colorscale[cnt]
                cnt += 1
                cumu_color_dict[idx] = trace_color
            trace_opacity = 1.0
            trace_info = {
                'name': str(idx), 
                'x': val[display_ndces['x']], 
                'y': val[display_ndces['y']], 
                'selectedpoints': spoint_ndces_this_trace, 
                'hoverinfo': 'text+name', 
                'text': point_ids_this_trace, 
                'mode': 'markers', 
                'opacity': trace_opacity, 
                'marker': {
                    'size': marker_size if str(idx) != 'Unknown' else bg_marker_size, 
                    'opacity': app_config.params['marker_opacity_factor'] if str(idx) != 'Unknown' else app_config.params['bg_marker_opacity_factor'], 
                    'symbol': 'circle', 
                    'color': trace_color
                }, 
                'selected': style_selected
            }
            if not app_config.params['three_dims']:
                trace_info.update({'type': 'scattergl'})
            else:
                trace_info.update({ 'type': 'scatter3d', 'z': val[display_ndces['z']] })
            traces_list.append(trace_info)
    return traces_list


def layout_scatter(annots):
    display_ndces = app_config.params['display_coordinates']
    new_layout = building_block_divs.create_scatter_layout(annots)
    return new_layout


def build_main_scatter(data_df, color_var, colorscale, highlight=False, 
                       bg_marker_size=app_config.params['bg_marker_size_factor'], 
                       marker_size=app_config.params['marker_size_factor'], 
                       annots=[], selected_point_ids=[], 
                       style_selected = building_block_divs.style_selected
                      ):
    if highlight:
        style_selected['marker']['color'] = 'white'
    else:
        style_selected['marker'].pop('color', None)    # Remove color if exists
    trace_list = traces_scatter(
        data_df, 
        color_var, 
        colorscale, 
        selected_point_ids, 
        bg_marker_size=bg_marker_size, 
        marker_size=marker_size, 
        style_selected=style_selected
    )
    return { 
        'data': trace_list, 
        'layout': layout_scatter(annots)
    }



# ========================================================
# =================== Raw data heatmap ===================
# ========================================================


def hm_row_scatter(fit_data, scatter_fig, hm_point_names, view_cocluster):
    row_scat_traces = []
    all_hm_point_names = []
    hmscat_mode = 'markers'
    # Decide if few enough points are around to display row labels
    if len(hm_point_names) <= 35:
        hmscat_mode = 'markers+text'
    if (scatter_fig is not None) and ('data' in scatter_fig):
        pts_so_far = 0
        # Re-sort rows or not? should be >=1 trace, so can check if first is continuous
        resort_rows = (not view_cocluster)
        is_continuous_color = None
        hm_row_ndces = []
        for trace in scatter_fig['data']:
            trace_markers = trace['marker']
            if 'line' in trace_markers:
                trace_markers['line']['width'] = 0.0
            if is_continuous_color is None:
                is_continuous_color = isinstance(trace_markers['color'], (list, tuple, np.ndarray))
            # Of the point names, choose the ones in this trace and get their indices...
            hm_point_names_this_trace = np.intersect1d(hm_point_names, trace['text'])           # point IDs in this trace
            num_in_trace = len(hm_point_names_this_trace)
            hm_point_ndces_this_trace = np.where(np.isin(hm_point_names, hm_point_names_this_trace))[0]        # this trace's row indices in heatmap
            y_coords_this_trace = np.arange(len(hm_point_names))[hm_point_ndces_this_trace]
            
            # At this point, rows are sorted in order of co-clustering. 
            # Now sort points by color within each trace. 
            # This does nothing if the colors are discrete (many traces), and is just for continuous plotting.
            if resort_rows:   # Row order determined by sorting, and continuous variable being plotted
                y_coords_this_trace = np.arange(pts_so_far, pts_so_far+num_in_trace)
                if is_continuous_color:
                    # Extract subset of rows in heatmap
                    resorted_ndces_this_trace = np.argsort(
                        np.array(trace_markers['color'])[hm_point_ndces_this_trace])
                    trace_markers['color'] = np.array(trace_markers['color'])[resorted_ndces_this_trace]
                    # TODO y_coords_this_trace = y_coords_this_trace[resorted_ndces_this_trace]
                else:
                    pts_so_far += num_in_trace
                    resorted_ndces_this_trace = np.arange(len(hm_point_names_this_trace))
                hm_point_names_this_trace = hm_point_names_this_trace[resorted_ndces_this_trace]
                hm_point_ndces_this_trace = hm_point_ndces_this_trace[resorted_ndces_this_trace]
            
            hm_row_ndces.extend(hm_point_ndces_this_trace)
            all_hm_point_names.extend(hm_point_names_this_trace)
            new_trace = {
                'name': trace['name'], 
                'x': np.zeros(num_in_trace), 
                'y': y_coords_this_trace, 
                'xaxis': 'x2', 
                'hoverinfo': 'text+name', 
                'text': hm_point_names_this_trace, 
                'mode': hmscat_mode, 
                'textposition': 'center left', 
                'textfont': building_block_divs.hm_font_macro, 
                'marker': trace_markers, 
                'selected': building_block_divs.style_selected, 
                'type': 'scatter'
            }
            row_scat_traces.append(new_trace) 
        # reorganize matrix if things were re-sorted.
        if resort_rows:
            fit_data = np.array([]) if len(hm_row_ndces) == 0 else fit_data[np.array(hm_row_ndces), :]
    return row_scat_traces, fit_data, all_hm_point_names


def hm_hovertext(data, rownames, colnames):
    pt_text = []
    # First the rows, then the cols
    for r in range(data.shape[0]):
        pt_text.append(["Gene: {}".format(str(rownames[r])) for k in data[r, :]])
        for c in range(data.shape[1]):
            pt_text[r][c] += "<br>Cell line: {}<br>Essentiality score: {}".format(str(colnames[c]), str(round(data[r][c], 3)))
    return pt_text


def display_heatmap_cb(
    hm_raw_data,    # 2D numpy array of selected data
    feat_names,     # col labels of hm_raw_data
    hm_point_names,    # (unique!) row labels of hm_raw_data
    scatter_fig,    # Scatterplot panel which this is mirroring.
    view_cocluster,  
    scatter_frac_domain=0.10
):
    fit_data = hm_raw_data
    if not app_config.params['hm_diverging']:
        fit_data = rna_log1pnorm(fit_data)
    # Identify (interesting) genes to plot. Currently: high-variance genes
    feat_ndces = interesting_feat_ndces(fit_data)
    absc_labels = feat_names[feat_ndces]
    fit_data = fit_data[:, feat_ndces]
    # Quantile normalize the data if necessary to better detect patterns.
    if app_config.params['hm_qnorm_plot']:
        qtiles = np.zeros_like(fit_data)
        nnz_ndces = np.nonzero(fit_data)
        qtiles[nnz_ndces] = sp.stats.rankdata(fit_data[nnz_ndces]) / len(fit_data[nnz_ndces])
        fit_data = qtiles
    # Spectral coclustering to cluster the heatmap. We always order rows (points) by spectral projection, 
    # But cols (features) can have different orderings for different viewing options.
    if (fit_data.shape[0] > 1):
        ordered_rows, ordered_cols = dm.compute_coclustering(fit_data)
        fit_data = fit_data[ordered_rows, :]
        hm_point_names = hm_point_names[ordered_rows]
    else:
        ordered_cols = np.arange(fit_data.shape[1])
    fit_data = fit_data[:, ordered_cols]
    absc_labels = absc_labels[ordered_cols]
    
    # Copy trace metadata from scatter_fig, in order of hm_point_names, to preserve colors etc.
    row_scat_traces, fit_data, hm_point_names = hm_row_scatter(fit_data, scatter_fig, hm_point_names, view_cocluster)
    pt_text = hm_hovertext(fit_data, hm_point_names, absc_labels)
    hm_trace = {
        'z': fit_data, 
        'x': absc_labels, 
        # 'y': hm_point_names, 
        'hoverinfo': 'text',
        'text': pt_text, 
        'colorscale': app_config.params['hm_colorscale'],
        'zmin': 0, 
        'colorbar': {
            'len': 0.3, 
            'thickness': 20, 
            'xanchor': 'left', 
            'yanchor': 'top', 
            'title': 'Ess. score',
            'titleside': 'top',
            'ticks': 'outside', 
            'titlefont': building_block_divs.colorbar_font_macro, 
            'tickfont': building_block_divs.colorbar_font_macro
        }, 
        'type': 'heatmap'
    }
    if app_config.params['hm_diverging']:
        max_magnitude = np.percentile(np.abs(fit_data), 99) if fit_data.shape[0] > 0 else 2
        hm_trace['zmin'] = -max_magnitude
        hm_trace['zmax'] = max_magnitude
    return {
        'data': [ hm_trace ] + row_scat_traces, 
        'layout': building_block_divs.create_hm_layout(scatter_frac_domain) 
    }


def generate_percluster_viz(raw_data, cell_cluster_list, cell_color_list, featID='Gene'):
    cluster_IDs = np.unique(cell_cluster_list)
    plot_values = np.zeros((len(cluster_IDs), raw_data.shape[1]))
    for i in range(len(cluster_IDs)):
        plot_values[i, :] = np.mean(raw_data[cell_cluster_list == cluster_IDs[i], :], axis=0)
    panel_layout = {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 0}, 
        'hovermode': 'closest', 
        'autosize': True,
        'xaxis': {
            'title': {'text': featID, 'font': building_block_divs.legend_font_macro }, 
            'tickfont': building_block_divs.legend_font_macro 
        }, 
        'yaxis': {
            'showticklabels': False,
            'automargin': True, 
            'ticks': 'outside', 
            'tickcolor': app_config.params['legend_font_color']
        },
        'plot_bgcolor': app_config.params['bg_color'],
        'paper_bgcolor': app_config.params['bg_color'], 
        'showlegend': True, 
        'legend': {
            'font': building_block_divs.legend_font_macro
        }
    }
    go_results = np.array(gp.gprofile(selected_genes))
    top_go_logpvals = np.array([])
    top_go_termnames = np.array([])
    top_go_dbIDs = np.array([])
    if go_results.shape[0] > 0:
        go_results = go_results[:topk, :]
        top_go_logpvals = -np.log10(go_results[:,2].astype(float))
        top_go_dbIDs = go_results[:,9]
        top_go_termnames = go_results[:,11]
    database_colors = { 'MF': '#CB3C19'}
    database_IDs = {'MF': 'Molecular function'}
    bar_colors = np.array([database_colors[x] for x in top_go_dbIDs])
    panel_data = []
    ordi = np.arange(len(top_go_dbIDs))[::-1]
    for c in np.unique(bar_colors):
        trace_locs = np.where(bar_colors == c)[0]
        trace_locs = trace_locs[::-1]      # Reverse because it's better.
        panel_data.append({
            'name': database_IDs[top_go_dbIDs[trace_locs[0]]], 
            'x': top_go_logpvals[trace_locs],
            'y': ordi[trace_locs], 
            'hovertext': [
                "-log(p): {}<br>{}<br>{}".format(
                    str(round(top_go_logpvals[t], 2)), 
                    top_go_termnames[t]
                ) for t in trace_locs], 
            'text': top_go_termnames[trace_locs], 
            'hoverinfo': 'text', 
            'insidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'outsidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'marker': { 'color': c },      
            'textposition': 'auto', 
            'orientation': 'h', 
            'type': 'bar'
        })
    return {'data': panel_data, 'layout': panel_layout }




from gprofiler import GProfiler
"""
Update GO enrichment panel.

g:GOSt API (in class header of gprofiler.py).

* ``all_results`` - (*Boolean*) All results, including those deemed not significant.
* ``ordered`` - (*Boolean*) Ordered query.
* ``exclude_iea`` - (*Boolean*) Exclude electronic GO annotations.
* ``underrep`` - (*Boolean*) Measure underrepresentation.
* ``evcodes`` - (*Boolean*) Request evidence codes in output as the
  final column.
* ``hier_sorting`` - (*Boolean*) Sort output into subgraphs.
* ``hier_filtering`` - (*Boolean*) Hierarchical filtering.
* ``max_p_value`` - (*Float*) Custom p-value threshold.
* ``min_set_size`` - (*Int*) Minimum size of functional category.
* ``max_set_size`` - (*Int*) Maximum size of functional category.
* ``min_isect_size`` - (*Int*) Minimum size of query / functional
  category intersection.
* ``max_isect_size`` - (*Int*) Maximum size of query / functional
  category intersection.
* ``correction_method`` - Algorithm used for multiple testing correction, one of:
  - ``GProfiler.THR_GSCS`` **Default** g:SCS.
  - ``GProfiler.THR_FDR`` Benjamini-Hochberg FDR.
  - ``GProfiler.THR_BONFERRONI`` Bonferroni.
* ``domain_size`` - Statistical domain size, one of:
  - ``GProfiler.DOMAIN_ANNOTATED`` - **Default** Only annotated genes.
  - ``GProfiler.DOMAIN_KNOWN`` - All known genes.
* ``custom_bg`` - (*String* | *List*) Custom statistical background
* ``src_filter`` - (*List*) A list of data source ID strings, e.g.
  ``["GO:BP", "KEGG"]``. These currently include GO (GO:BP, GO:MF,
  GO:CC to select a particular GO branch), KEGG, REAC, TF, MI, CORUM,
  HP, HPA, OMIM.
"""
def display_goenrich_panel_func(selected_genes, topk=20):
    panel_layout = {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 25}, 
        'hovermode': 'closest',
        'orientation': 90,
        'autosize': True,
        'xaxis': {
            'title': {'text': '-log(p)', 'font': building_block_divs.legend_font_macro }, 
            'tickfont': building_block_divs.legend_font_macro 
        }, 
        'yaxis': {
            'showticklabels': False,
            'automargin': True, 
            'ticks': 'outside', 
            'tickcolor': app_config.params['legend_font_color']
        },
        'plot_bgcolor': app_config.params['bg_color'],
        'paper_bgcolor': app_config.params['bg_color'], 
        'showlegend': True, 
        'legend': {
            'font': building_block_divs.legend_font_macro
        }
    }
    gp = GProfiler("MyToolName/0.1")
    go_results = np.array(gp.gprofile(selected_genes))
    top_go_logpvals = np.array([])
    top_go_termnames = np.array([])
    top_go_dbIDs = np.array([])
    if go_results.shape[0] > 0:
        go_results = go_results[:topk, :]
        top_go_logpvals = -np.log10(go_results[:,2].astype(float))
        top_go_dbIDs = go_results[:,9]
        top_go_termnames = go_results[:,11]
        top_go_termIDs = go_results[:,8]
        top_go_queryhits = go_results[:,13]
    database_colors = { 'MF': '#CB3C19', 'BP': '#FA9D18', 'CC': '#198520', 'keg': '#D9869B', 'rea': '#335ACC', 'tf': '#4F6E9B', 'mir': '#53B8AD', 'hpa': '#542DB1', 'cor': '#6AAB19', 'hp': '#8C0784'}
    database_IDs = {'MF': 'Molecular function', 'BP': 'Biological process', 'CC': 'Cellular component', 'keg': 'KEGG', 'rea': 'REAC', 'tf': 'TF', 'mir': 'MIRNA', 'hpa': 'HPA', 'cor': 'CORUM', 'hp': 'HP'}
    bar_colors = np.array([database_colors[x] for x in top_go_dbIDs])
    panel_data = []
    ordi = np.arange(len(top_go_dbIDs))[::-1]
    for c in np.unique(bar_colors):
        trace_locs = np.where(bar_colors == c)[0]
        trace_locs = trace_locs[::-1]      # Reverse because it's better.
        panel_data.append({
            'name': database_IDs[top_go_dbIDs[trace_locs[0]]], 
            'x': top_go_logpvals[trace_locs],
            'y': ordi[trace_locs], 
            # 'y': top_go_termIDs[trace_locs],
            'hovertext': [
                "-log(p): {}<br>{}<br>{}".format(
                    str(round(top_go_logpvals[t], 2)), 
                    top_go_termnames[t], 
                    top_go_termIDs[t]
                ) for t in trace_locs], 
            'text': top_go_termnames[trace_locs], 
            'hoverinfo': 'text', 
            'insidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'outsidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'marker': { 'color': c },      
            'textposition': 'auto', 
            'orientation': 'h', 
            'type': 'bar'
        })
    return {'data': panel_data, 'layout': panel_layout }
