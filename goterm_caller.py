''' Python port of g:ProfileR, the R wrapper for the g:Profiler functional enrichment tool.
From https://github.com/vals/python-gprofiler
'''

import requests
import pandas as pd
from goatools.associations import read_ncbi_gene2go
from goatools.go_search import GoSearch



BASE_URL = "http://biit.cs.ut.ee/gprofiler/"
HEADERS = {'User-Agent': 'python-gprofiler'}

def gprofiler(query, organism='hsapiens', ordered_query=False, significant=True,
              exclude_iea=False, region_query=False, max_p_value=1.0, max_set_size=0,
              correction_method='analytical', hier_filtering='none',
              domain_size='annotated', custom_bg=[], numeric_ns='', no_isects=False,
              png_fn=None, include_graph=False, src_filter=None, 
              mode='enrich'
             ):
    ''' Annotate gene list functionally
    Interface to the g:Profiler tool for finding enrichments in gene lists. Organism
    names are constructed by concatenating the first letter of the name and the
    family name. Example: human - 'hsapiens', mouse - 'mmusculus'. If requesting PNG
    output, the request is directed to the g:GOSt tool in case 'query' is a vector
    and the g:Cocoa (compact view of multiple queries) tool in case 'query' is a
    list. PNG output can fail (return FALSE) in case the input query is too large.
    In such case, it is advisable to fall back to a non-image request.
    Returns a pandas DataFrame with enrichment results or None
    '''
    query_url = ''
    
    mode = 'enrich'   # Currently gprofiler gconvert gives an HTTP 405 error.
    if mode == 'enrich':
        my_url = BASE_URL + 'gcocoa.cgi'
    elif mode == 'lookup':
        my_url = BASE_URL + 'convert.cgi'
    wantpng = True if png_fn else False
    output_type = 'mini_png' if wantpng else 'mini'

    if wantpng:
        raise NotImplementedError('PNG Output not implemented')
        return

    if include_graph:
        raise NotImplementedError('Biogrid Interactions not implemented (include_graph)')
        return

    # Query

    qnames = list(query)

    if not qnames:
        raise ValueError('Missing query')

    query_url = ' '.join(qnames)

    # Significance threshold

    if correction_method == 'gSCS':
        correction_method = 'analytical'

    if correction_method not in ('analytical', 'fdr', 'bonferroni'):
        raise ValueError("Multiple testing correction method not recognized (correction_method)")

    # Hierarchical filtering

    if hier_filtering not in ('none', 'moderate', 'strong'):
        raise ValueError("hier_filtering must be one of \"none\", \"moderate\" or \"strong\"")

    if hier_filtering == 'strong':
        hier_filtering = 'compact_ccomp'
    elif hier_filtering == 'moderate':
        hier_filtering = 'compact_rgroups'
    else:
        hier_filtering = ''

    # Domain size

    if domain_size not in ('annotated', 'known'):
        raise ValueError("domain_size must be one of \"annotated\" or \"known\"")

    # Custom background

    if isinstance(custom_bg, list):
        custom_bg = ' '.join(custom_bg)
    else:
        raise TypeError('custom_bg need to be a list')

    # Max. set size

    if max_set_size < 0:
        max_set_size = 0

    # HTTP request

    if mode == 'enrich':
        query_params = {
            'organism': organism,
            'query': query_url,
            'output': output_type,
            'analytical': '1',
            'sort_by_structure': '1',
            'ordered_query': '1' if ordered_query else '0',
            'significant': '1' if significant else '0',
            'no_iea': '1' if exclude_iea else '0',
            'as_ranges': '1' if region_query else '0',
            'omit_metadata': '0' if include_graph else '1',
            'user_thr': str(max_p_value),
            'max_set_size': str(max_set_size),
            'threshold_algo': correction_method,
            'hierfiltering': hier_filtering,
            'domain_size_type': domain_size,
            'custbg_file': '',
            'custbg': custom_bg,
            'prefix': numeric_ns,
            'no_isects': '1' if no_isects else '0'
        }
    elif mode == 'lookup':
        query_params = {
            'query': query_url, 
            'target': 'GO'
        }

    if src_filter:
        for i in src_filter:
            query_params['sf_' + i] = '1'


    raw_query = requests.post(my_url, data=query_params, headers=HEADERS)

    # Here PNG request parsing would go, but not implementing that

    if wantpng:
        pass

    # Requested text

    split_query = raw_query.text.split('\n')

    # Here interaction parsing would go, but not implementing that

    if include_graph:
        pass

    # Parse main result body

    split_query = [
        s.split('\t')
        for s in split_query
        if s and not s.startswith('#')
    ]
        
    enrichment = pd.DataFrame(split_query)
    if mode == 'enrich':
        colnames = [
            "query.number", "significant", "p.value", "term.size", 
            "query.size", "overlap.size", "recall", "precision", 
            "term.id", "domain", "subgraph.number", "term.name", 
            "relative.depth", "intersection"
        ]
        numeric_colnames = [
            "query.number", "p.value", "term.size", 
            "query.size", "overlap.size", "recall", 
            "precision", "subgraph.number", "relative.depth"
        ]
    elif mode == 'lookup':
        return enrichment
        colnames = [
            "alias.number", "alias", "target.number", 
            "target", "name", "description", "namespace"
        ]
        numeric_colnames = [
            "alias.number", "target.number"
        ]
    
    if enrichment.shape[1] > 0:
        enrichment.columns = colnames
        enrichment.index = enrichment['term.id']
        numeric_columns = numeric_colnames
        for column in numeric_columns:
            enrichment[column] = pd.to_numeric(enrichment[column])
        if mode == 'enrich':
            enrichment['significant'] = enrichment['significant'] == '!'
    else:
        enrichment = None

    return enrichment



"""
The following is an adaptation of https://github.com/tanghaibao/goatools/blob/master/goatools/go_search.py
"""

import sys
from goatools.obo_parser import GODag

class GoSearch(object):
    """Returns GOs matching a regex pattern."""

    def __init__(self, fin_go_basic_obo, go2items, log=None, verbose=False):
        self.verbose = verbose
        self.log = sys.stdout if log is None else log
        self.bstdout = True if log is None else log
        # Some obo fields often used in searching. Many are optional to load when reading obo
        self.goa_srch_hdrs = ['defn', 'comment', 'name', 'is_a', 'relationship', 'synonym', 'xref']
        self.obo_dag = GODag(fin_go_basic_obo, optional_attrs=self.goa_srch_hdrs)
        self.go2items = go2items

    def get_matching_gos(self, compiled_pattern, **kws):
        """Return all GOs which match the user regex pattern."""
        # kws: prt gos
        matching_gos = []
        obo_dag = self.obo_dag
        prt = kws['prt'] if 'prt' in kws else self.log
        if self.verbose:
            prt.write('\nPATTERN SEARCH: "{P}"\n'.format(P=compiled_pattern.pattern))
        # Only look through GOs in annotation or user-specified GOs
        srchgos = kws['gos'] if 'gos' in kws else self.go2items.keys()
        for go_id in srchgos:
            go_obj = obo_dag.get(go_id, None)
            if go_obj is not None:
                for hdr in self.goa_srch_hdrs:
                    if hdr in go_obj.__dict__:
                        fld_val = getattr(go_obj, hdr)
                        matches = self._search_vals(compiled_pattern, fld_val)
                        for mtch in matches:
                            if self.verbose:
                                prt.write("MATCH {go_id}({NAME}) {FLD}: {M}\n".format(
                                    FLD=hdr, go_id=go_obj.id, NAME=go_obj.name, M=mtch))
                        if matches:
                            matching_gos.append(go_id)
            else:
                prt.write("**WARNING: {GO} found in annotation is not found in obo\n".format(
                    GO=go_id))
        matching_gos = set(matching_gos)
        # Print summary message
        if self.verbose:
            self._summary_matching_gos(prt, compiled_pattern.pattern, matching_gos, srchgos)
        return matching_gos

    @staticmethod
    def _summary_matching_gos(prt, pattern, matching_gos, all_gos):
        """Print summary for get_matching_gos."""
        msg = 'Found {N} GO(s) out of {M} matching pattern("{P}")\n'
        num_gos = len(matching_gos)
        num_all = len(all_gos)
        prt.write(msg.format(N=num_gos, M=num_all, P=pattern))

    def _search_vals(self, compiled_pattern, fld_val):
        """Search for user-regex in scalar or iterable data values."""
        matches = []
        if isinstance(fld_val, set):
            for val in fld_val:
                self._search_val(matches, compiled_pattern, val)
        elif isinstance(fld_val, str):
            self._search_val(matches, compiled_pattern, fld_val)
        return matches

    @staticmethod
    def _search_val(matches, compiled_pattern, fld_val):
        """Search for user-regex in scalar data values."""
        mtch = compiled_pattern.search(fld_val)
        if mtch:
            matches.append(fld_val)

    def add_children_gos(self, gos):
        """Return children of input gos plus input gos."""
        lst = []
        obo_dag = self.obo_dag
        get_children = lambda go_obj: list(go_obj.get_all_children()) + [go_obj.id]
        for go_id in gos:
            go_obj = obo_dag[go_id]
            lst.extend(get_children(go_obj))
        return set(lst)

    def get_items(self, gos):
        """Given GO terms, return genes or gene products for the GOs."""
        items = []
        for go_id in gos:
            items.extend(self.go2items.get(go_id, []))
        return set(items)