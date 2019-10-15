    # -*- coding: utf-8 -*-
import numpy as np, pandas as pd

params = {}

params['title'] = "Gene co-essentiality browser"


_DEPLOY_LOCALLY = False

if not _DEPLOY_LOCALLY:
    data_pfx = '/mnt/lab_data/kundaje/abalsubr/coessentiality/'
    # data_pfx = '/srv/www_coessentiality/coessentiality-browser/'
else:
    data_pfx = '/Users/akshay/github/coessentiality-browser/'

param_mixes = [
    "{}data/vizdf{}.csv".format(data_pfx, sffix) for sffix in ["_GLS01_CO99", "_GLS1_CO0", "_GLS02_CO49_RC49"]
]

clusterone_params = ['0.2', '0.5', '0.8', '0.9', '0.95', '0.99', '0.8_batch_corrected', '0.9_batch_corrected', '0.95_batch_corrected', '0.99_batch_corrected']
params['plot_data_df_path'] = param_mixes + [data_pfx + "data/vizdf_GLS01_CO99_d_{}.csv".format(x) for x in clusterone_params]
params['raw_ess_data_path'] = "{}data/essentiality.tsv.gz".format(data_pfx)
params['mutation_data_path'] = "{}data/CCLE_DepMap_18q3_maf_20180718.txt".format(data_pfx)
params['mutation_arr_path'] = "{}data/CCLE_DepMap_18q3_mutations.npy".format(data_pfx)
params['gene_fullnames_path'] = "{}data/hgnc_gene_symbol_to_name.tsv".format(data_pfx)
params['genenames_path'] = "{}data/gene_long_names.tsv".format(data_pfx)
params['gene_ensID_path'] = "{}data/gene_ensIDs.tsv".format(data_pfx)
params['shRNA_data_path'] = "{}data/D2_combined_gene_dep_scores.csv".format(data_pfx)
params['expression_data_path'] = "{}data/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct".format(data_pfx)
params['expression_arr_path'] = "{}data/CCLE_DepMap_18q3_RNAseq_RPKM.npy".format(data_pfx)
params['expression_cell_lines_path'] = "{}data/CCLE_DepMap_18q3_cell_lines.npy".format(data_pfx)
params['video_path'] = "{}data/45sec_browser_video.mp4".format(data_pfx)

params['string_ppi_names_path'] = "{}data/9606.protein.info.v11.0.txt.gz".format(data_pfx)
params['string_ppi_network_path'] = "{}data/9606.protein.links.v11.0.txt.gz".format(data_pfx)
params['string_ppi_matrix_path'] = "{}data/string_ppi_mat.npz".format(data_pfx)
params['string_matrix_ascoess_path'] = "{}data/string_mat.npz".format(data_pfx)

params['cheng_matrix_path'] = "{}data/cheng_mat.npz".format(data_pfx)

# params['dataset_options'] = [x.split('/')[-1].split('.')[0] for x in params['plot_data_df_path']]
podata = ["ClusterONE clusters (d = {})".format(x) for x in clusterone_params]
params['dataset_options'] = ["ClusterONE clusters", "GLS", "Mixed"] + podata

params['gene2go_path'] = data_pfx + 'gene2go'
params['go_obo_path'] = data_pfx + 'data/go-basic.obo'
params['gotermnames_path'] = data_pfx + 'gotermnames.npy'
params['gotermIDs_path'] = data_pfx + 'gotermIDs.npy'

params['bg_color'] = '#000000'

params['hm_colorvar_name'] = 'Ess.'
params['hm_diverging'] = True
params['default_color_var'] = 'Colors'



params['display_coordinates'] = { 'x': 'hUMAP_x',  'y': 'hUMAP_y' }
params['qnorm_plot'] = False
params['hm_qnorm_plot'] = False
params['continuous_color'] = False




# ======================================================
# ================== Colors and sizes ==================
# ======================================================

# Custom colorscales.

cmap_celltypes = ["#f7ff00","#fabebe","#ff8300","#f000ff","#74ee15","#4363d8","#001eff","#ff3300","#cab2d6","#008080","#808000","#bcf60c","#bec1d4","#fb9a99","#ffd8b1","#3cb44b","#bc5300","#ffe119","#33ccff","#911eb4","#46f0f0","#d220c8","#e6beff","#e6194b","#aaffc3","#000075"]

# Perceptually uniform blackbody, good for black background as in https://github.com/kennethmoreland-com/kennethmoreland-com.github.io/blob/master/color-advice/black-body/black-body.ipynb
cmap_custom_blackbody = [[0.0, "#000000"], [0.39, "#b22222"], [0.58, "#e36905"], [0.84, "#eed214"], [1.0, "#ffffff"]]

# From https://github.com/BIDS/colormap/blob/master/parula.py
# pc = [matplotlib.colors.to_hex(x) for x in parulac]; d = np.arange(len(pc)); d = np.round(d/max(d), 4); parula = [x for x in zip(d, pc)]
cmap_parula = [(0.0, '#352a87'), (0.0159, '#363093'), (0.0317, '#3637a0'), (0.0476, '#353dad'), (0.0635, '#3243ba'), (0.0794, '#2c4ac7'), (0.0952, '#2053d4'), (0.1111, '#0f5cdd'), (0.127, '#0363e1'), (0.1429, '#0268e1'), (0.1587, '#046de0'), (0.1746, '#0871de'), (0.1905, '#0d75dc'), (0.2063, '#1079da'), (0.2222, '#127dd8'), (0.2381, '#1481d6'), (0.254, '#1485d4'), (0.2698, '#1389d3'), (0.2857, '#108ed2'), (0.3016, '#0c93d2'), (0.3175, '#0998d1'), (0.3333, '#079ccf'), (0.3492, '#06a0cd'), (0.3651, '#06a4ca'), (0.381, '#06a7c6'), (0.3968, '#07a9c2'), (0.4127, '#0aacbe'), (0.4286, '#0faeb9'), (0.4444, '#15b1b4'), (0.4603, '#1db3af'), (0.4762, '#25b5a9'), (0.4921, '#2eb7a4'), (0.5079, '#38b99e'), (0.5238, '#42bb98'), (0.5397, '#4dbc92'), (0.5556, '#59bd8c'), (0.5714, '#65be86'), (0.5873, '#71bf80'), (0.6032, '#7cbf7b'), (0.619, '#87bf77'), (0.6349, '#92bf73'), (0.6508, '#9cbf6f'), (0.6667, '#a5be6b'), (0.6825, '#aebe67'), (0.6984, '#b7bd64'), (0.7143, '#c0bc60'), (0.7302, '#c8bc5d'), (0.746, '#d1bb59'), (0.7619, '#d9ba56'), (0.7778, '#e1b952'), (0.7937, '#e9b94e'), (0.8095, '#f1b94a'), (0.8254, '#f8bb44'), (0.8413, '#fdbe3d'), (0.8571, '#ffc337'), (0.873, '#fec832'), (0.8889, '#fcce2e'), (0.9048, '#fad32a'), (0.9206, '#f7d826'), (0.9365, '#f5de21'), (0.9524, '#f5e41d'), (0.9683, '#f5eb18'), (0.9841, '#f6f313'), (1.0, '#f9fb0e')]

# Default discrete colormap for <= 20 categories, from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/. See also http://phrogz.net/css/distinct-colors.html and http://tools.medialab.sciences-po.fr/iwanthue/
cmap_custom_discrete = ["#bdbdbd", 
                        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']#, '#7d87b9', '#bec1d4', '#d6bcc0']
# Custom red/blue diverging for black background, from https://gka.github.io/palettes
cmap_custom_rdbu_diverging = [[0.0, '#0000ff'], [0.1111, '#442dfa'], [0.2222, '#6b59e0'], [0.3333, '#6766a3'], [0.4444, '#323841'], [0.5555, '#483434'], 
                              [0.6666, '#b3635b'], [0.7777, '#ee5d49'], [0.8888, '#ff3621'], [1.0, '#ff0000']]
# Custom yellow/blue diverging for black background. From the following code:
# x = sns.diverging_palette(227, 86, s=98, l=77, n=20, center='dark').as_hex(); [s for s in zip(np.arange(len(x))/(len(x)-1), x)]
cmap_custom_ylbu_diverging = [(0.0, '#3acdfe'), (0.05263157894736842, '#37bbe6'), (0.10526315789473684, '#35a9cf'), (0.15789473684210525, '#3295b6'), (0.21052631578947367, '#2f829e'), (0.2631578947368421, '#2d6f85'), (0.3157894736842105, '#2a5d6e'), (0.3684210526315789, '#274954'), (0.42105263157894735, '#25373d'), (0.47368421052631576, '#222324'), (0.5263157894736842, '#232322'), (0.5789473684210527, '#363621'), (0.631578947368421, '#474720'), (0.6842105263157895, '#5a5a1e'), (0.7368421052631579, '#6b6b1d'), (0.7894736842105263, '#7e7e1c'), (0.8421052631578947, '#8f901b'), (0.8947368421052632, '#a2a21a'), (0.9473684210526315, '#b3b318'), (1.0, '#c4c417')]
cmap_custom_orpu_diverging = [(0.0, '#c2b5fe'), (0.05263157894736842, '#b1a5e6'), (0.10526315789473684, '#a096cf'), (0.15789473684210525, '#8e85b6'), (0.21052631578947367, '#7c759e'), (0.2631578947368421, '#6a6485'), (0.3157894736842105, '#59556e'), (0.3684210526315789, '#464354'), (0.42105263157894735, '#35343d'), (0.47368421052631576, '#232324'), (0.5263157894736842, '#242323'), (0.5789473684210527, '#3d332a'), (0.631578947368421, '#544132'), (0.6842105263157895, '#6e523a'), (0.7368421052631579, '#856041'), (0.7894736842105263, '#9e7049'), (0.8421052631578947, '#b67f50'), (0.8947368421052632, '#cf8f58'), (0.9473684210526315, '#e79d5f'), (1.0, '#feac66')]


if 'colorscale_discrete' not in params:
    params['colorscale_discrete'] = cmap_custom_discrete
if 'colorscale_continuous' not in params:
    params['colorscale_continuous'] = cmap_custom_ylbu_diverging    # 'Viridis'

if params['continuous_color']:
    params['colorscale'] = params['colorscale_continuous']
else:
    params['colorscale'] = params['colorscale_discrete']

if params['hm_diverging']:
    params['hm_colorscale'] = cmap_custom_ylbu_diverging
else:
    params['hm_colorscale'] = cmap_custom_blackbody


params['hover_edges'] = ""
params['edge_color'] = 'rgb(255,255,255)'
params['edge_width'] = 1
params['incl_edges'] = False
params['three_dims'] = False

params['bg_marker_size_factor'] = 3.6
params['marker_size_factor'] = 3.5
params['bg_marker_opacity_factor'] = 0.5
params['marker_opacity_factor'] = 1.0
params['font_color'] = 'white'

params['legend_bgcolor'] = '#000000'
params['legend_bordercolor'] = 'white'
# params['legend_borderwidth'] = 2
params['legend_font_color'] = 'white'
params['legend_font_size'] = 16
params['hm_font_size'] = 6

params['upload_img_path'] = data_pfx + 'assets/upload.png'
params['download_img_path'] = data_pfx + 'assets/download.png'

# Some things are best left to depend on the size of the data - opacity changes with number of points plotted!
