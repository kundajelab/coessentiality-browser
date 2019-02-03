"""
Routines for efficiently computing eigenmaps on graphs. 
Author: Akshay Balsubramani
"""
import numpy as np, os, time, scipy as sp, sklearn
from sklearn import preprocessing
from sklearn.cluster.bicluster import SpectralCoclustering
from scipy.sparse import issparse

"""
'TerminalInteractiveShell': terminal IPython
'ZMQInteractiveShell': Jupyter (notebook AND qtconsole)
Fails (NameError): regular Python interpreter
s/o https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook/39662359#39662359
"""
def isnotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


def compute_coclustering(
    fit_data, 
    num_clusters=1, 
    tol_bicluster=0.005,  # sparsity otherwise annoyingly causes underflows w/ sklearn
):
    if num_clusters == 1:
        num_clusters = min(fit_data.shape[0], 5)
    model = SpectralCoclustering(n_clusters=num_clusters, random_state=0)
    model.fit(fit_data + tol_bicluster)
    ordered_rows = np.argsort(model.row_labels_)
    ordered_cols = np.argsort(model.column_labels_)
    return ordered_rows, ordered_cols


"""
An instance of this object stores the information about a particular alignment of the data. 
Typically this is a coordinate system obtained by analyzing variation of the combined graph. 
Input: two graph adjacency matrices.
"""
class sc_aligner(object):
    
    def __init__(self, graph_list, mode='pls'):
        #self.common_dims = common_dims
        self.graph_list = graph_list
        self.corresp_pairs = []
        if mode == 'pls':
            pass
        else:
            self.build_combined_graph()
    
    # Build the graph with optional correspondence info
    def build_combined_graph(self, corresp=[], corresp_reg=[]):
        num_cells_total = np.sum([g.shape[0] for g in graph_list])
        self.full_graph = sp.sparse.csr_matrix(np.zeros((num_cells_total, num_cells_total)))
        dimsofar = 0
        # The block diagonal terms...
        for g in self.graph_list:
            thisdim = g.shape[0]
            self.full_graph[(dimsofar):(dimsofar+thisdim), (dimsofar):(dimsofar+thisdim)] = g
            dimsofar += thisdim
        # The cross (correspondence) terms...
        # TODO: Generalize this to multiway correspondences by adding a weight to the graph.
        for i in range(len(corresp)):
            subset1 = corresp[i][0]
            subset2 = corresp[i][1]     # Both assumed to be numpy arrays of length >= 1
            self.full_graph[subset1, subset2] = g
            dimsofar += thisdim
        return
    
    # Calculates (eigen)decomposition of matrix as required.
    def decompose(self, train_pos_ndces=None, train_neg_ndces=None, random_state=None):
        print("Finding top {} eigenvectors...".format(common_dims))
        train_mats = []



# ===============================================================
# =============== Diffusion eigenmap computations ===============
# ===============================================================

# def is_symmetric(A):
#     return np.allclose(A, A.T)

"""
A couple of helper functions for decomposing a given pairwise similarity matrix.
"""
def compute_transitions(adj_mat, alpha=1.0, sym=True):
    """Compute symmetrized transition matrix. 
    alpha : The exponent of the diffusion kernel. 
    * 1 = Laplace-Beltrami (density)-normalized [default]. 
    * 0.5 = normalized graph Laplacian (Fokker-Planck dynamics, cf. "Diffusion maps, spectral clustering and eigenfunctions of Fokker-Planck operators" NIPS2006).
    * 0 = classical graph Laplacian. 
    """
    similarity_mat = symmetric_part(adj_mat)
    dens = np.asarray(similarity_mat.sum(axis=0))  # dens[i] is an estimate for the sampling density at point i.
    K = sp.sparse.spdiags(np.power(dens, -alpha), 0, similarity_mat.shape[0], similarity_mat.shape[0])
    W = sp.sparse.csr_matrix(K.dot(similarity_mat).dot(K))
    z = np.sqrt(np.asarray(W.sum(axis=0)).astype(np.float64))    # sqrt(density)
    Zmat = sp.sparse.spdiags(1.0/z, 0, W.shape[0], W.shape[0])
    return Zmat.dot(W).dot(Zmat) if sym else Zmat.power(2).dot(W)


def compute_eigen(adj_mat, n_comps=2, sym=True):
    """
    Compute eigendecomposition of sparse transition kernel matrix.
    (Computing in 64-bit can matter for analysis use beyond just visualization!).
    sym indicates that the transition matrix should be symmetrized.
    NOTE: use np.linalg.eigh if we want to support non-sparse matrices.
    """
    if sym:
        eigvals, evecs = sp.sparse.linalg.eigsh(adj_mat.astype(np.float64), k=n_comps)
    else:
        # eigvals, evecs = sp.sparse.linalg.eigs(adj_mat.astype(np.float64), k=n_comps)   # DON'T USE without further thresholding: complex-number underflow issues
        evecs, eigvals, _ = sp.sparse.linalg.svds(adj_mat.astype(np.float64), k=n_comps)
    #eigvals, evecs = eigvals.astype(np.float32), evecs.astype(np.float32)
    sorted_ndces = np.argsort(np.abs(eigvals))[::-1]
    return eigvals[sorted_ndces], evecs[:, sorted_ndces]


def heat_eigval_dist(eigvals, t):
    return np.power(np.abs(eigvals), t)/np.sum(np.power(np.abs(eigvals), t))


""" Returns the symmetrized version of the input matrix. (The anti-symmetrized version is A - symmetric_part(A) .) """
def symmetric_part(A):
    return 0.5*(A + A.T)


"""
Compute diffusion map embedding.
sym_compute: Whether to compute the SVD on a symmetric matrix.
sym_return: Whether to return the SVD of the symmetrized transition matrix.
Assumes that {sym_return => sym_compute} holds, i.e. doesn't allow sym_compute==False and sym_return==True.
Returns (n_comps-1) components, excluding the first which is the same for all data.
"""
def diffmap_proj(adj_mat, 
                 t=None, 
                 min_energy_frac=0.95, 
                 n_dims=None,      # Number of dims to return
                 n_comps=None,     # Number of comps to use in computation.
                 return_eigvals=False, 
                 embed_type='diffmap', 
                 sym_compute=True, 
                 sym_return=False):
    if n_comps is None:     # When data are high-d dimensional with log2(d) \leq 14-16 as for scRNA, 2K eigenvector computation is tolerably fast; change otherwise
        n_comps = min(2000, adj_mat.shape[0]-1)
    if n_dims is None:
        n_dims = n_comps - 1
    eigvals, evecs = compute_eigen(compute_transitions(adj_mat, sym=sym_compute), n_comps=n_comps, sym=sym_compute)
    if sym_compute:
        evecs_sym = evecs
        evecs_unsym = np.multiply(evecs, np.outer(1.0/evecs[:,0].astype(np.float64), np.ones(evecs.shape[1])))
    else:
        evecs_unsym = evecs
    if sym_return:
        if not sym_compute:
            print("TODO: ERROR LOGGED HERE")
            return
        eigvecs_normalized = preprocessing.normalize(np.real(evecs_sym), axis=0, norm='l2')
    else:
        eigvecs_normalized = preprocessing.normalize(np.real(evecs_unsym), axis=0, norm='l2')
    
    if t is None:     # Use min_energy_frac to determine the fraction of noise variance 
        t = min_t_for_energy(eigvals, n_dims+1, min_energy_frac)
    frac_energy_explained = np.cumsum(np.power(np.abs(eigvals), t)/np.sum(np.power(np.abs(eigvals), t)))[n_dims]
    print("{} dimensions contain about {} fraction of the variance in the first {} dimensions (Diffusion time = {})".format(
        n_dims+1, frac_energy_explained, n_comps, t))
    if embed_type=='naive':
        all_comps = eigvecs_normalized
    elif embed_type=='diffmap':
        all_comps = np.power(np.abs(eigvals), t) * eigvecs_normalized
    elif embed_type=='commute': # Check the correctness!
        all_comps = np.power((1-np.abs(eigvals)), -t/2) * eigvecs_normalized
    if not return_eigvals:
        return all_comps[:,1:(n_dims+1)]     # Return n_dims dimensions, skipping the first trivial one.
    else:
        return (all_comps[:,1:(n_dims+1)], eigvals)    # Return the eigenvalues as well.


def min_t_for_energy(eigvals, desired_dim, min_energy_frac, max_t=None):
    # Calc upper bound for t principal eigengap (if g=lbda1/lbda2, then g^t < 100 implies t < log(100)/log(g) ). Don't change unless you know what you're doing!
    if max_t is None:
        max_t = np.log(100)/np.log(max(eigvals[0]/eigvals[1], 1.01))
    f = lambda t: (np.sum(heat_eigval_dist(eigvals, t)[:desired_dim]) - min_energy_frac)
    if f(0)*f(max_t) >= 0:    # since f is always nondecreasing this means the zero isn't in the interval
        return max_t if f(0) < 0 else 0
    return sp.optimize.brentq(f, 0, max_t)


"""
* SLOW, NAIVE IMPLEMENTATION - WHEN UMAP ISN'T AVAILABLE
Build kNN graph, returning a sparse matrix. First, we build the matrix accretively from zero, cycling through the nodes in some arbitrary order and setting the kNN's of each to be nonzero. So the degree of each vertex is \geq k. Call the matrix A.
Then symmetrize (or not) A in some way:
- 'mutual': min(A, A^T), making degree \leq k
- 'inclusive': max(A, A^T), making degree \geq k
- 'matching': TODO solve an OT problem! Called b-matching in the literature.
See \cite{JebaraWC09}. 
"""
def build_knn(mat, k=10, symmetrize_type='inclusive'):
    sparse_adj = sp.sparse.csr_matrix(np.zeros(mat.shape))
    for i in range(mat.shape[0]):
        matarr = mat[i,:].toarray().flatten()
        nbrs = np.argsort(matarr)[::-1][:k]    # Highest k similarities
        sparse_adj[i, nbrs] = 1.0#matarr[nbrs]
    if (symmetrize_type == 'mutual'):
        return sparse_adj.minimum(sparse_adj.transpose())
    elif (symmetrize_type == 'inclusive'):
        return sparse_adj.maximum(sparse_adj.transpose())
    else:
        print("Mode not yet implemented.")
        return sparse_adj


"""
Compute local geometry descriptors in terms of the rate of local heat dissipation. 
Input: a set of points (rows) in their diffusion-component representation (possibly just the main few components), with eigenvalues.
Returns: Local geometry descriptors, one for each point.
"""
def local_descriptor(reduced_data, eigvals, filtertype='HKS'):
    sqnorms = np.multiply(reduced_data, reduced_data)
    if (filtertype == 'HKS'):
        return sqnorms
    return




# =================================================================
# =============== Learning alignment transformation ===============
# =================================================================





# ======================================================================
# =============== Graph operations (modified from PyGSP) ===============
# ======================================================================
# If unfamiliar with this area, see https://arxiv.org/pdf/1408.5781 for a useful primer.

# Spectral coclustering of the rows and columns of a dense numpy matrix.
# This can be interpreted as the min-cut solution of a bipar
def spectral_coclustering(adj_mat):
    return


# Return incidence matrix of graph as list. Each edge between i and j gets counted twice, from i->j and j->i. 
# The outgoing vertex is marked positive, and the incoming negative.
# Input: Sparse adjacency matrix.
def graph_as_list(adj_mat):
    v_in, v_out = adj_mat.nonzero()
    weights = np.asarray(adj_mat[v_in, v_out]).squeeze()
    return v_in, v_out, weights


# Calculate given pairwise similarity kernel over a given subset of edges.
def calc_over_edges(graph_in, data_in_values, dist_fn=None):
    gg = np.zeros_like(graph_in.toarray())
    nnzlist = list(zip(*np.nonzero(graph_in.toarray())))
    for i in range(len(nnzlist)):
        d = nnzlist[i]
        if dist_fn is None:
            gg[d] = np.exp(-np.sum(np.square(data_in_values[d[0],:] - data_in_values[d[1],:])))
        else:
            gg[d] = dist_fn(data_in_values[d[0],:], data_in_values[d[1],:])
    return sp.sparse.csr_matrix(gg)


# Return the sparse (edges x nodes) directed incidence matrix - the graph differential operator.
def del_op(adj_mat):
    v_in, v_out, weights = graph_as_list(adj_mat)
    num_edges = len(v_in)
    rows = np.concatenate((np.arange(num_edges), np.arange(num_edges)))
    cols = np.empty(2*num_edges)
    cols[:num_edges] = v_in
    cols[num_edges:] = v_out
    vals = np.empty(2*num_edges)
    vals[:num_edges] = np.sqrt(weights)
    vals[num_edges:] = -vals[:num_edges]
    return sp.sparse.csc_matrix((vals, (rows, cols)), shape=(num_edges, adj_mat.shape[0]))


# Divergence of a function in a neighborhood on the data manifold. Operates on a vector field (function over edges)
def div_op_graph(graph_arg, fn):
    return del_op(graph_arg).transpose().dot(fn)


# Gradient of a function over the graph vertices.
def grad_op_graph(graph_arg, fn):
    return del_op(graph_arg).dot(fn)


# Laplacian of a function over the graph.
def laplacian_op_graph(graph_arg, fn):
    return del_op(graph_arg).transpose().dot(del_op(graph_arg)).dot(fn)


# Return smoothed version of function (by rolling forward the random walk a bit.)
# Can smooth several signals at once by stacking them together as columns of a matrix.
def local_smooth_op(graph_arg, fn, frac_stationary=0.5):
    T = compute_transitions(graph_arg, sym=False)
    nbr_vals = T.dot(fn)
    return (frac_stationary*fn + (1-frac_stationary)*nbr_vals)


# Computes curvature for each edge, returning matrix of such curvatures.
def ricci_curvature(graph_arg, alpha=0.8):
    # Compute (1-alpha)-diffused distributions for each vertex.
    orig_sig = np.identity(graph_arg.shape[0])
    smoodist = local_smooth_op(graph_arg, orig_sig, frac_stationary=alpha)
    # TODO: compute EMD between each pair of distributions.
    allpairs_emd = None
    # Compute all-pairs graph distances and thereby Ricci curvature.
    allpairs_dist = sp.sparse.csgraph.floyd_warshall(graph_arg, directed=True)
    return 1 - np.divide(allpairs_emd, allpairs_dist)
    


# Modulate spectrum according to a function.
# TODO include different times.
def wavelet_op(graph_arg, base_fn, t=1, init_dist=None):
    newtrans = matrix_lift(base_fn, graph_arg)
    if init_dist is None:
        return newtrans
    else:
        return newtrans.dot(init_dist)


# Compute trend filtering denoised version of fn over graph.
def trendfilter_smooth(graph_arg, fn):
    return


# Apply fn pointwise to matrix's spectrum 
def matrix_lift(fn, matrix):
    # Find eigendecomposition.
    # Apply fn to eigvals pointwise.
    return  # re-aggregated matrix


"""
Returns MLE estimator of local intrinsic dimension at each data point.
Poisson process axioms: (1) nonoverlapping regions are independent, (2) Prob of event is proportional to interval size, (3) <=1 event occurs in an interval.
Closed under thinning, mapping, and displacement (by transition matrix). Fit inhomogeneous too.
TODO: What we really want is a measure of canalization, i.e. robustness to natural trajectory variation. If we can predict the effect of gene products in this canalized state, we can uncover the TFs underpinning the state.
"""
def intrinsic_dim_MLE(weighted_kNN_graph, k, impute_1nn=False, bias_correction=True, eps=1e-5):
    # TODO: Use entropy of this dist. with k to detect levels of structure. Add other estimation methods if necessary.
    # weighted_kNN_graph is assumed to be similarities = e^{-distances} where distances are provided.
    idim_fn = np.zeros(weighted_kNN_graph.shape[0])
    if issparse(weighted_kNN_graph):
        g = weighted_kNN_graph.toarray()
    else:
        g = weighted_kNN_graph
    wtknn_dists = []
    idim_fn = np.zeros(g.shape[0])
    for i in range(g.shape[0]):
        knn_dists = np.sort(g[i,:])[-1:-(k+1):-1]
        r = np.log(-(np.log(knn_dists) + eps))
        if impute_1nn:    # If the similarity to the most similar point is already 1, pretend it's as similar as the 2nd most
            r[0] = r[1]
        idim_fn[i] = 1.0/np.mean(r[-1] - r)
    if bias_correction:
        idim_fn = idim_fn*((k-2)/(k-1))
    return idim_fn


# Return a (target_dim x source_dim) linear projection matrix U where U^T U = I, i.e. the cols are an orthonormal basis.
# X U is a low-dimensional representation of X that preserves pairwise L2-norms. If projections == True, return this projected data instead of basis.
# TODO add diffmap embedding
# TODO use target_energy to auto-choose components
def dimreduce_basis(target_dim, source_dim, data=None, target_energy=0.0, mode='random', projections=False, random_state=None):
    if mode == 'random':
        toret = np.random.random(size=(source_dim, target_dim))
        basis = sp.linalg.orth(toret)
    elif mode == 'random_sparse':
        basis = sklearn.random_projection.sparse_random_matrix(n_components=target_dim, n_features=source_dim, density='auto', random_state=random_state).transpose()
    elif mode == 'diff_map':
        return "Not yet implemented"     # TODO proper centralized error logging
    elif mode == 'pca':
        # Assume sparse count data, so don't zero-center!
        mpc = sc.pp.pca(data, zero_center=False, n_comps=target_dim, return_info=True)
        basis = mpc[1].transpose()
    return basis.transpose().dot(data.transpose()).transpose() if projections else basis


# Pairwise computation of mean FPTs. Warning: overly dependent on local geometry!
def calc_mean_first_passage_time(transition_mat, stationary_dist):
    dim = transition_mat.shape[0]
    W_mat = np.outer(np.ones(dim), stationary_dist)
    fund_mat = scipy.mat.inv(np.eye() - transition_mat + 1)
    fpt = 0 # Should be n_jj - n_ij/(w_j) where n is fund_mat
    return fpt




# =======================================================
# =============== Visualization utilities ===============
# =======================================================


"""
import matplotlib.pyplot as plt

# Plot first dim eigenvalues of two spectra against each other, scaling time appropriately so they contain the same energy.

def plot_spectra(eigvals1, eigvals2, dim, min_energy_frac=0.95):
    t1 = min_t_for_energy(eigvals1, dim, min_energy_frac=min_energy_frac)
    t2 = min_t_for_energy(eigvals2, dim, min_energy_frac=min_energy_frac)
    absc = np.power(np.abs(eigvals1)[:dim], t1)
    ordi = np.power(np.abs(eigvals2)[:dim], t2)
    plt.scatter(absc, ordi, marker='x')
    #plt.scatter(-np.log(absc), -np.log(ordi) )
    plt.xlabel('Batch 1 eigenvalues')
    plt.ylabel('Batch 2 eigenvalues')
    plt.plot(np.arange(0.0,1,0.01), np.arange(0.0,1,0.01), linestyle='--', color='r')
    plt.show()
    return


# Plot histogram of k^{th} diffusion component projection of each point. 
def plot_histogram(unnorm_projections, k=1, numbins=150):
    dmat = unnorm_projections * np.sqrt(unnorm_projections.shape[0])    # Each axis projection now has zero mean, unit variance.
    plt.hist(dmat[:,k-1], bins=numbins, range=(-2,2))
    plt.show()
    return
"""



# =================================================================
# =============== Metrics for goodness of alignment ===============
# =================================================================

"""
Statistic for determining goodness of alignment of dataset A to B; a two-way kBET.
"""
def kbet_2batch(adj_mat_A, adj_mat_B):
    return


"""
Calculate and read off the distance to the k^th nearest neighbor, as a vector over nodes. Assumes each point has zero distance to self.
"""
def calc_knn_dists(adj_mat, k=5):
    nbrdists = np.zeros(adj_mat.shape[0])
    for i in range(adj_mat.shape[0]):
        nbrdists[i] = np.sort(adj_mat[i,:])[k-1]
    return nbrdists


"""
# TODO: Finish by calling https://github.com/HazyResearch/hyperbolics
def hyperbolic_embedding():
    return

# TODO: Landscape viz (potential plotted on dim-reduced plane of single points). Use https://moderndata.plot.ly/3d-surface-plots-in-ipython-notebook-and-plotly/
def waddfig4viz():
    return
"""


"""
If DataStructure() is a:
Stack - DFS
Queue - BFS

TODO: Fix this by specifying tree (node) API.
"""
def x_first_search(start):
    seen = {start}
    pending = DataStructure()
    pending.push(start)
    while pending:
        current = pending.pop()
        for neighbor in current.neighbors:
            if neighbor not in seen:
                seen.add(neighbor)
                pending.push(neighbor)



# ============================================
# =============== Misc metrics ===============
# ============================================

# Construct gene-gene graph from cluster memberships; edges are Jaccard similarities 
def pairwise_jaccard_graph(input_memberships):
    a = input_memberships.dot(input_memberships.T).toarray()
    tmp = np.add(np.add(-a.T, np.diagonal(a)).T, np.diagonal(a))
    return sp.sparse.csr_matrix(np.nan_to_num(np.divide(a, tmp)))


# Tukey fences


# Rolling statistics in numpy, quickly:
# http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html



"""
Extend kernel to other points in the space besides the given data. Drop some dimensions if possible, i.e. if their energy doesn't bottleneck estimation of the function on a held-out subset. Jensen-Shannon between ground truth and heldout at those locations.
"""
def extend_nystrom():
    return
