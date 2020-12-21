# Application-specific routines for working with (smallish matrices of) data.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd, diffmap as dm
import app_config, building_block_divs, goterm_caller
# from goatools.associations import read_ncbi_gene2go
# from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT
import re
import umap
from sklearn.decomposition import TruncatedSVD



# =======================================================
# ================== Utility functions ==================
# =======================================================

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



# =========================================================
# =================== Main scatter plot ===================
# =========================================================

"""
(Data, layout) for the main graph panel.
Color_var is either a field of the plotting df, or a numpy array.
"""


full_gene_names = np.array(pd.read_csv(app_config.params['genenames_path'], sep="\t", header=None, index_col=False)).flatten()


# Here selected_point_ids is a list of unique string IDs of points. 
def traces_scatter(
    data_df, 
    color_var, 
    colorscale, 
    selected_point_ids, 
    bg_marker_size=app_config.params['bg_marker_size_factor'], 
    marker_size=app_config.params['marker_size_factor'], 
    style_selected=building_block_divs.style_selected, 
    point_names=None
):
    traces_list = []
    display_ndces = app_config.params['display_coordinates']
    cumu_color_dict = {}
    # Check to see if color_var is continuous or discrete and plot points accordingly
    if isinstance(color_var, (list, tuple, np.ndarray)):     # Color_var is an array, not a col index.
        continuous_color_var = color_var
        spoints = np.where(np.isin(point_names, selected_point_ids))[0]
        colorbar_title = app_config.params['hm_colorvar_name']
        pt_text = ["{}<br>{}<br>Essentiality score: {}".format(
            point_names[i], full_gene_names[i], round(continuous_color_var[i], 3)) for i in range(len(point_names))]
        max_magnitude = np.percentile(np.abs(continuous_color_var), 99)
        traces_list.append({ 
            'name': 'Data', 
            'x': data_df[display_ndces['x']], 
            'y': data_df[display_ndces['y']], 
            'selectedpoints': spoints, 
            'hoverinfo': 'text', 
            'hovertext': pt_text, 
            'text': point_names, 
            'mode': 'markers', 
            'marker': {
                'size': bg_marker_size, 
                'opacity': app_config.params['marker_opacity_factor'], 
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
                'colorscale': app_config.cmap_custom_ylbu_diverging, 
                'cmin': -max_magnitude, 
                'cmax': max_magnitude
            }, 
            'selected': style_selected, 
            'type': 'scattergl'
        })
    else:    # Categorical color scheme, one trace per color
        cnt = 0
        for idx, val in data_df.groupby(color_var):
            point_ids_this_trace = list(val['gene_names'])
            spoint_ndces_this_trace = np.where(np.isin(point_ids_this_trace, selected_point_ids))[0]
            if idx not in cumu_color_dict:
                trace_color = colorscale[cnt]
                cnt += 1
                cumu_color_dict[idx] = trace_color
            trace_opacity = 1.0
            full_ids_this_trace = full_gene_names#[np.isin(point_names, point_ids_this_trace)]
            pt_text = ["{}<br>{}".format(point_ids_this_trace[i], full_ids_this_trace[i]) for i in range(len(point_ids_this_trace))]
            # pt_text = ["{}".format(point_ids_this_trace[i]) for i in range(len(point_ids_this_trace))]
            trace_info = {
                'name': str(idx), 
                'x': val[display_ndces['x']], 
                'y': val[display_ndces['y']], 
                'selectedpoints': spoint_ndces_this_trace, 
                'hoverinfo': 'text', 
                'hovertext': pt_text, 
                'text': point_ids_this_trace, 
                'mode': 'markers', 
                'opacity': trace_opacity, 
                'marker': {
                    'size': marker_size if str(idx) in ['Unknown', 'Unannotated'] else bg_marker_size, 
                    'opacity': app_config.params['marker_opacity_factor'], 
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


def build_main_scatter(
    data_df, color_var, colorscale, highlight=False, 
    bg_marker_size=app_config.params['bg_marker_size_factor'], 
    marker_size=app_config.params['marker_size_factor'], 
    annots=[], selected_point_ids=[], 
    style_selected = building_block_divs.style_selected, 
    point_names=None
):
    if highlight:
        style_selected['marker']['color'] = '#ff4f00' # Golden gate bridge red
        style_selected['marker']['size'] = 10     # TODO: Change this back on unhighlighting.
    else:
        style_selected['marker'].pop('color', None)    # Remove color if exists
    trace_list = traces_scatter(
        data_df, 
        color_var, 
        colorscale, 
        selected_point_ids, 
        bg_marker_size=bg_marker_size, 
        marker_size=marker_size, 
        style_selected=style_selected, 
        point_names=point_names
    )
    display_ndces = app_config.params['display_coordinates']
    return { 
        'data': trace_list, 
        'layout': building_block_divs.create_scatter_layout(annots)
    }



# ========================================================
# =================== Raw data heatmap ===================
# ========================================================


def hm_row_scatter(scatter_fig, hm_point_names, view_cocluster, row_clustIDs=None):
    row_scat_traces = []
    all_hm_point_names = []
    hmscat_mode = 'markers'
    # Decide if few enough points are around to display row labels
    if len(hm_point_names) <= 30:
        hmscat_mode = 'markers+text'
    if (scatter_fig is not None) and ('data' in scatter_fig):
        pts_so_far = 0
        is_continuous_color = None
        for trace in scatter_fig['data']:
            trace_markers = trace['marker']
            trace_markers['showscale'] = False    # Masks colorbar
            if 'line' in trace_markers:
                trace_markers['line']['width'] = 0.0
            if is_continuous_color is None:
                is_continuous_color = isinstance(trace_markers['color'], (list, tuple, np.ndarray))
            # Of the point names, choose the ones in this trace and get their indices...
            hm_point_where_this_trace = np.isin(hm_point_names, trace['text'])
            hm_point_names_this_trace = hm_point_names[hm_point_where_this_trace]
            num_in_trace = len(hm_point_names_this_trace)
            hm_point_ndces_this_trace = np.where(hm_point_where_this_trace)[0]        # this trace's row indices in heatmap
            y_coords_this_trace = np.arange(len(hm_point_names))[hm_point_ndces_this_trace]
            # At this point, rows are sorted in order of co-clustering. 
            all_hm_point_names.extend(hm_point_names_this_trace)
            new_trace = {
                'name': trace['name'], 
                'x': np.zeros(num_in_trace), 
                'y': y_coords_this_trace, 
                'xaxis': 'x2', 
                'hoverinfo': 'text', 
                'text': hm_point_names_this_trace, 
                'mode': hmscat_mode, 
                'textposition': 'center left', 
                'textfont': building_block_divs.hm_font_macro, 
                'marker': trace_markers, 
                'selected': building_block_divs.style_selected, 
                'type': 'scatter'
            }
            row_scat_traces.append(new_trace)
    return row_scat_traces, all_hm_point_names


def hm_col_plot(
    fit_data, 
    feat_colordict, 
    reordered_groups=None, 
    reordered_featnames=None, 
    col_clustIDs=None, 
    geneview_mode='Mutation', 
    geneview_gene=None, 
    geneview_data=None, 
    geneview_celllines=None
):
    if reordered_featnames is None:
        reordered_featnames = feat_colordict.keys()
    if reordered_groups is None:
        reordered_groups = reordered_featnames
    transform_col_order = np.arange(fit_data.shape[1])
    fit_data = fit_data[:, transform_col_order]
    reordered_groups = reordered_groups[transform_col_order]
    reordered_featnames = reordered_featnames[transform_col_order]
    # Plot the gene viewer, integrating Depmap mutation/expression databases.
    col_scat_traces = []
    if geneview_gene is not None:
        pt_text = []
        gvdata = np.zeros_like(reordered_featnames)    # np.zeros((1, len(reordered_featnames)))
        if geneview_mode == 'Mutation':
            mutdata = geneview_data[geneview_data[:, 0] == geneview_gene, :]
            mutdata_celllines = mutdata[:, 6]
            where_mutations = np.isin(reordered_featnames, mutdata_celllines)
            for i in range(len(reordered_featnames)):
                newtext = "Cell line: {}".format(reordered_featnames[i])
                if where_mutations[i]:
                    query_results = mutdata[(mutdata_celllines == reordered_featnames[i]), :]
                    gvdata[i] = query_results.shape[0]
                    qwe = ["chr{}:{}    {}    Protein: {}    Codon: {}".format(
                        query_results[i, 1], query_results[i, 2], query_results[i, 4], 
                        query_results[i, 11], query_results[i, 10]#, query_results[i, 7], query_results[i, 8], query_results[i, 9]
                    ) for i in range(query_results.shape[0])]
                    newtext = newtext + "<br>{}".format("<br>".join(qwe))
                pt_text.append(newtext)
            gvdata = np.minimum(gvdata, 1)
            to_extend = [{
                'y': gvdata, 'x': reordered_featnames, 'yaxis': 'y3', 
                'hoverinfo': 'text', 'hovertext': pt_text, 'hoverdistance': 40, 
                'colorscale': app_config.cmap_custom_blackbody, 
                'marker': { 'color': 'white' }, 'type': 'bar'
            }]
        elif geneview_mode == 'Expression':
            exprdata_celllines = geneview_celllines
            exprdata_genes = geneview_data[:, 0]
            # where_exprdata = np.isin(reordered_featnames, exprdata_celllines)
            print(np.mean(where_exprdata), len(reordered_featnames), len(exprdata_celllines))
            for i in range(len(reordered_featnames)):
                newtext = "Cell line: {}".format(reordered_featnames[i])
                if reordered_featnames[i] in exprdata_celllines:
                    gvdata[i] = geneview_data[(exprdata_genes == geneview_gene), (exprdata_celllines == reordered_featnames[i])][0]
                    newtext = newtext + "<br>Expression: {}".format(gvdata[i])
                pt_text.append(newtext)
            to_extend = [{
                'y': np.ones(len(gvdata)), 
                'x': reordered_featnames, #'y': np.reshape(gvdata, (1, len(gvdata))), 
                'yaxis': 'y3', 
                'hoverinfo': 'text', 'hovertext': pt_text, 'hoverdistance': 40, 
                'colorbar': {
                    'len': 0.3, 
                    'thickness': 20, 
                    'xanchor': 'left', 
                    'yanchor': 'top', 
                    'title': 'Expression',
                    'titleside': 'top',
                    'ticks': 'outside', 
                    'titlefont': building_block_divs.colorbar_font_macro, 
                    'tickfont': building_block_divs.colorbar_font_macro
                }, 
                'marker': { 'colorscale': app_config.cmap_custom_blackbody, 'color': list(np.log(np.array(gvdata).astype(float) + 1)) }, 
                'type': 'bar'
            }]
        col_scat_traces.extend(to_extend)
    # Make scatterplot of cell lines.
    for feat_group in feat_colordict:
        trace_marker = { 'color': feat_colordict[feat_group] }
        feat_ndces_this_trace = np.where(reordered_groups == feat_group)[0]
        feat_names_this_trace = reordered_featnames[feat_ndces_this_trace]
        x_coords_this_trace = feat_names_this_trace# np.arange(len(reordered_featnames))[feat_ndces_this_trace]
        new_trace = {
            'name': feat_group, 
            'x': x_coords_this_trace, 
            'y': np.zeros(len(feat_names_this_trace)), 
            'yaxis': 'y2', 
            'hoverinfo': 'text+name', 
            'text': feat_names_this_trace, 
            'mode': 'markers', 
            'textposition': 'top center', 
            'textfont': building_block_divs.hm_font_macro, 
            'marker': trace_marker, 
            'selected': building_block_divs.style_selected, 
            'type': 'scatter'
        }
        col_scat_traces.append(new_trace)
    return col_scat_traces, fit_data


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
    feat_colordict={}, 
    feat_group_names=None, 
    geneview_mode='Mutation', 
    geneview_gene=None, 
    geneview_data=None, 
    geneview_celllines=None, 
    scatter_frac_domain=0.13, 
    scatter_frac_range=0.08, 
    show_legend=False
):
    fit_data = hm_raw_data
    if fit_data is None or len(fit_data.shape) < 2:
        return
    # Identify (interesting) cell lines to plot. Currently: high-variance ones
    feat_ndces = interesting_feat_ndces(fit_data)
    absc_labels = feat_names[feat_ndces]
    absc_group_labels = feat_group_names[feat_ndces]
    fit_data = fit_data[:, feat_ndces]
    
    # Quantile normalize the data if necessary to better detect patterns.
    if app_config.params['hm_qnorm_plot']:
        qtiles = np.zeros_like(fit_data)
        nnz_ndces = np.nonzero(fit_data)
        qtiles[nnz_ndces] = sp.stats.rankdata(fit_data[nnz_ndces]) / len(fit_data[nnz_ndces])
        fit_data = qtiles
    # Spectral coclustering to cluster the heatmap. We always order rows (points) by spectral projection, but cols (features) can have different orderings for different viewing options.
    row_clustIDs = np.zeros(fit_data.shape[0])
    col_clustIDs = np.zeros(fit_data.shape[1])
    if (fit_data.shape[0] > 1):
        ordered_rows, ordered_cols, row_clustIDs, col_clustIDs = dm.compute_coclustering(fit_data)
        fit_data = fit_data[ordered_rows, :]
        hm_point_names = hm_point_names[ordered_rows]
    else:
        ordered_cols = np.arange(fit_data.shape[1])
    fit_data = fit_data[:, ordered_cols]
    absc_labels = absc_labels[ordered_cols]
    
    if absc_group_labels is not None:
        absc_group_labels = absc_group_labels[ordered_cols]
    # Copy trace metadata from scatter_fig, in order of hm_point_names, to preserve colors etc.
    row_scat_traces = []
    row_scat_traces, hm_point_names = hm_row_scatter(scatter_fig, hm_point_names, view_cocluster, row_clustIDs=row_clustIDs)
    col_scat_traces, fit_data = hm_col_plot(
        fit_data, 
        feat_colordict, 
        reordered_groups=absc_group_labels, 
        reordered_featnames=absc_labels, 
        col_clustIDs=col_clustIDs, 
        geneview_mode=geneview_mode, 
        geneview_gene=geneview_gene, 
        geneview_data=geneview_data, 
        geneview_celllines=geneview_celllines
    )
    pt_text = hm_hovertext(fit_data, hm_point_names, absc_labels)
    hm_trace = {
        'z': fit_data, 
        'x': absc_labels, 
        'customdata': hm_point_names, 
        'hoverinfo': 'text',
        'text': pt_text, 
        'colorscale': app_config.params['hm_colorscale'], 
        'colorbar': {
            'len': 0.3, 
            'thickness': 20, 
            'xanchor': 'left', 
            'yanchor': 'top', 
            'title': 'Less essential',
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
    # assemble coordinates of lines adumbrating clusters.
    clustersep_line_coords = []
    for cid in np.unique(col_clustIDs):
        ndcesc = np.where(col_clustIDs == cid)[0]
        clustersep_line_coords.append(np.min(ndcesc) - 0.5)
    if len(fit_data.shape) > 1:
        clustersep_line_coords.append(fit_data.shape[1] - 0.5)
    
    return {
        'data': [ hm_trace ] + row_scat_traces + col_scat_traces, 
        'layout': building_block_divs.create_hm_layout(
            scatter_frac_domain=scatter_frac_domain, 
            scatter_frac_range=scatter_frac_range, 
            show_legend=show_legend, 
            clustersep_coords=clustersep_line_coords, 
            yaxis_label=(len(hm_point_names) > 30)
        )
    }


# Given a gene set, returns GO+other database enrichment results using gProfiler.
def get_goenrichment_from_genes(gene_list):
    return goterm_caller.gprofiler(gene_list)


with open(app_config.data_pfx + "go2gene_dict.txt", "r") as f:
    w2 = f.read().strip()
go2gene_dict = json.loads(w2)

"""
# Given a GO term query, returns a combined list of genes uxnder that ID.
def get_genes_from_goterm(goterm_re_str, mode='gaf'):
    if len(goterm_re_str) == 0:
        return []
    if mode == 'regex':      # Given a regex, return using GO's association files; else given GO term ID(s).
        go2geneids_human = read_ncbi_gene2go(app_config.params['gene2go_path'], taxids=[9606], go2geneids=True)
        srchhelp = goterm_caller.GoSearch(app_config.params['go_obo_path'], go2items=go2geneids_human)
        gos = srchhelp.get_matching_gos(re.compile(goterm_re_str))
    else:
        gos = [goterm_re_str]
    return go2gene_dict[gos[0]] if gos[0] in go2gene_dict else []
"""


"""
Update GO enrichment panel.
https://biit.cs.ut.ee/gprofiler/page/apis. or 
g:GOSt API (in class header of gprofiler.py):
* ``all_results`` - (*Boolean*) All results, including those deemed not significant.
* ``ordered`` - (*Boolean*) Ordered query.
* ``exclude_iea`` - (*Boolean*) Exclude electronic GO annotations.
* ``underrep`` - (*Boolean*) Measure underrepresentation.
* ``evcodes`` - (*Boolean*) Request evidence codes in output as the final column.
* ``hier_sorting`` - (*Boolean*) Sort output into subgraphs.
* ``hier_filtering`` - (*Boolean*) Hierarchical filtering.
* ``max_p_value`` - (*Float*) Custom p-value threshold.
* ``min_set_size`` - (*Int*) Minimum size of functional category.
* ``max_set_size`` - (*Int*) Maximum size of functional category.
* ``min_isect_size`` - (*Int*) Minimum size of query / functional category intersection.
* ``max_isect_size`` - (*Int*) Maximum size of query / functional category intersection.
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
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'title': {'text': '-log(p)', 'font': building_block_divs.legend_font_macro }, 
            'tickfont': building_block_divs.legend_font_macro 
        }, 
        'yaxis': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
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
    if len(selected_genes) == 0:
        return {'data': [], 'layout': panel_layout }
    go_results = get_goenrichment_from_genes(selected_genes)
    top_go_logpvals = np.array([])
    top_go_termnames = np.array([])
    top_go_dbIDs = np.array([])
    if (go_results is not None) and (go_results.shape[0] > 0):
        x = np.array(np.argsort(go_results['p.value']))
        go_results = go_results.iloc[x[:topk], :]
        top_go_logpvals = -np.log10(go_results['p.value'].astype(float))
        top_go_dbIDs = go_results['domain']
        top_go_termnames = go_results['term.name']
        top_go_termIDs = go_results['term.id']
        top_go_queryhits = go_results['intersection']
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
            # 'y': ordi[trace_locs], 
            'y': top_go_termIDs[trace_locs],
            'hovertext': [
                "-log10(p): {}<br>{}<br>{}".format(
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
