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
    return (
        ordered_rows, ordered_cols, 
        model.row_labels_[ordered_rows], model.column_labels_[ordered_cols]
    )


def umap_embed_inplace(
    ann_heat, n_components = 2, min_dist = 0.5, spread = 1.0, n_epochs = 0, init_coords = 'spectral', 
    alpha = 1.0, gamma = 1.0, negative_sample_rate = 5, random_state = 0
):
    random_state = check_random_state(random_state)
    a, b = find_ab_params(spread, min_dist)
    X_umap = simplicial_set_embedding(
        ann_heat.X,
        ann_heat.uns['neighbors']['connectivities'].tocoo(),
        n_components, alpha, a, b, gamma, negative_sample_rate, n_epochs,
        init_coords, random_state, 'euclidean', {}, verbose=0
    )


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
    elif (symmetrize_type == 'fuzzy'):
        return sparse_adj.maximum(sparse_adj.transpose())
    else:
        print("Mode not yet implemented.")
        return sparse_adj


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


# ============================================
# =============== Misc metrics ===============
# ============================================

# Construct gene-gene graph from cluster memberships; edges are Jaccard similarities 
def pairwise_jaccard_graph(input_memberships):
    a = input_memberships.dot(input_memberships.T).toarray()
    tmp = np.add(np.add(-a.T, np.diagonal(a)).T, np.diagonal(a))
    return sp.sparse.csr_matrix(np.nan_to_num(np.divide(a, tmp)))

