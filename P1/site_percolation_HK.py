import numpy as np


# ============================================================
# Set USE_NUMBA = True to compile functions using Numba
# ============================================================

#USE_NUMBA = False ### Compila el codi. Per fer proves, no cal. Quan feu les simualcions llargues, poseu a True.
USE_NUMBA = True ### Compila el codi. Per fer proves, no cal. Quan feu les simualcions llargues, poseu a True.

# ============================================================
# Conditional Decorator for NUMBA 
# ============================================================

if (USE_NUMBA):
  try:
     import numba as nb
     njit = nb.njit
  except Exception:
     USE_NUMBA = False
     def njit(*args, **kwargs):
        # case: @njit without parentheses
        if args and callable(args[0]):
            return args[0]

        # case: @njit(...)
        def decorator(func):
            return func
        return decorator

else:
    def njit(*args, **kwargs):
        # case: @njit without parentheses
        if args and callable(args[0]):
            return args[0]

        # case: @njit(...)
        def decorator(func):
            return func
        return decorator


# ============================================================
# Numba-compileable Union-Find + build routine
# ============================================================

@njit(cache=True)
def uf_find(parent, x):
    # Path halving (iterative)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

@nb.njit(cache=True)
def uf_union(parent, size, a, b):
    ra = uf_find(parent, a)
    rb = uf_find(parent, b)
    if ra == rb:
        return
    if size[ra] < size[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    size[ra] += size[rb]

@nb.njit(cache=True)
def build_percolation(occ, L):
    """
    Parameters
    ----------
    occ : 1D boolean array of length N=L*L (flattened lattice)
    L   : int

    Returns
    -------
    sites        : 1D int64 array of occupied indices (flat)
    root_of_site : 1D int64 array length N; root id at occupied, -1 at empty
    counts       : 1D int64 array length N; cluster size per root id
    spans        : 1D boolean array length N; spans[root] True/False
    parent       : 1D int64 DSU parent (kept for completeness/debug)
    size         : 1D int64 DSU size (kept for completeness/debug)
    """
    N = L * L

    # --- DSU init ---
    parent = np.empty(N, dtype=np.int64)
    size = np.ones(N, dtype=np.int64)
    for i in range(N):
        parent[i] = i

    # --- occupied sites list ---
    # Numba supports np.nonzero / np.where
    sites = np.nonzero(occ)[0].astype(np.int64)

    # --- unions (right + down with PBC) ---
    for idx in range(sites.size):
        i = sites[idx]
        x = i // L
        y = i - x * L

        # right neighbor with PBC
        j = i + 1 if y != (L - 1) else i - (L - 1)
        if occ[j]:
            uf_union(parent, size, i, j)

        # down neighbor with PBC
        j = i + L if x != (L - 1) else i - L * (L - 1)
        if occ[j]:
            uf_union(parent, size, i, j)

    # --- root labeling for occupied sites + counts ---
    root_of_site = np.full(N, -1, dtype=np.int64)
    counts = np.zeros(N, dtype=np.int64)

    for idx in range(sites.size):
        i = sites[idx]
        r = uf_find(parent, i)
        root_of_site[i] = r
        counts[r] += 1

    # --- shadow spanning: unique roots per row / per column ---
    # Use "stamping" arrays to avoid np.unique.
    row_hits = np.zeros(N, dtype=np.int64)
    col_hits = np.zeros(N, dtype=np.int64)
    seen_row = np.full(N, -1, dtype=np.int64)  # seen_row[root] = last row stamp
    seen_col = np.full(N, -1, dtype=np.int64)  # seen_col[root] = last col stamp

    for idx in range(sites.size):
        i = sites[idx]
        r = root_of_site[i]  # root id
        x = i // L
        y = i - x * L

        if seen_row[r] != x:
            seen_row[r] = x
            row_hits[r] += 1

        if seen_col[r] != y:
            seen_col[r] = y
            col_hits[r] += 1

    spans = (row_hits == L) | (col_hits == L)

    return sites, root_of_site, counts, spans, parent, size


# ============================================================
# Public Python class (w API)
# ============================================================

class SitePercolation2D:
    """
    Numba-accelerated version of HK code.

    - PBC unions (right/down with wrap) for connectivity
    - "shadow spanning": a cluster spans if it has >=1 site in every row OR every column
    - cluster_matrix returns root IDs (0 for empty, root id for occupied) 
    - percolates returns (bool, size, root_id) where root_id matches cluster_matrix coding
    """

    def __init__(self, L: int, p: float, seed=None):

        self.L = int(L)
        self.p = float(p)
        self.rng = np.random.default_rng(seed)
        self.lattice = (self.rng.random((self.L, self.L)) < self.p)

        # cached build results
        self._built = False
        self._occ_sites = None
        self._root_of_site = None
        self._counts = None
        self._spans = None

        # optional DSU internals (useful if you want to inspect)
        self._parent = None
        self._size = None

    def rebuild(self, p=None, seed=None, lattice=None):
        """
        Rebuild with a new lattice (common in Monte Carlo runs).
        You can pass:
          - lattice: explicit boolean (L,L)
          - or p (and optional seed) to regenerate.
        """
        if lattice is not None:
            lattice = np.asarray(lattice, dtype=bool)
            if lattice.shape != (self.L, self.L):
                raise ValueError("lattice shape mismatch")
            self.lattice = lattice
        else:
            if p is not None:
                self.p = float(p)
            if seed is not None:
                self.rng = np.random.default_rng(seed)
            self.lattice = (self.rng.random((self.L, self.L)) < self.p)

        self._built = False
        self._occ_sites = self._root_of_site = self._counts = self._spans = None
        self._parent = self._size = None

    def _build(self):
        if self._built:
            return
        occ = self.lattice.ravel()
        sites, root_of_site, counts, spans, parent, size = build_percolation(occ, self.L)

        self._occ_sites = sites
        self._root_of_site = root_of_site
        self._counts = counts
        self._spans = spans
        self._parent = parent
        self._size = size
        self._built = True

    # --------------------------------------------------------
    # Percolation detection (shadow spanning)
    # --------------------------------------------------------
    def percolates(self):
        """
        Returns
        -------
        percolates : bool
        size       : int (size of spanning cluster; else max cluster size)
        root       : int (root id of spanning cluster; else root id of max)

        NOTE: 'root' matches the coding in cluster_matrix() (which outputs root ids).
        """
        self._build()

        counts = self._counts
        spans = self._spans

        existing = np.flatnonzero(counts)
        if existing.size == 0:
            return False, 0, -1

        # max cluster (fallback)
        rmax = int(existing[np.argmax(counts[existing])])
        smax = int(counts[rmax])

        span_roots = existing[spans[existing]]
        if span_roots.size:
            r = int(span_roots[0])
            return True, int(counts[r]), r

        return False, smax, rmax

    # --------------------------------------------------------
    # Cluster sizes list (called often)
    # --------------------------------------------------------
    def size_list(self, mode="all"):
        """
        mode:
          'all'         all clusters
          'perc'/'span' spanning clusters (shadow spanning)
          'fin'         non-spanning clusters
          'max'         size of maximum cluster only
        """
        self._build()

        counts = self._counts
        spans = self._spans

        existing = np.flatnonzero(counts)
        if existing.size == 0:
            return np.array([], dtype=np.int64)

        sizes = counts[existing]

        if mode == "all":
            return sizes.astype(np.int64, copy=True)

        if mode.startswith(("per", "spa", "inf")):
            return counts[existing[spans[existing]]].astype(np.int64, copy=True)

        if mode.startswith(("fin", "non", "sma")):
            return counts[existing[~spans[existing]]].astype(np.int64, copy=True)

        if mode.startswith(("max","tot")):
            return np.array([int(sizes.max())], dtype=np.int64)

        raise ValueError(f"Unknown mode: {mode!r}")

    # --------------------------------------------------------
    # Occasional labeled output
    # --------------------------------------------------------
    def cluster_matrix(self):
        """
        Returns an LxL int matrix:
          0 for empty
          root_id for occupied 
        """
        self._build()
        L = self.L
        labels_flat = np.zeros(L * L, dtype=np.int64)
        labels_flat[self._occ_sites] = self._root_of_site[self._occ_sites]
        return labels_flat.reshape((L, L))

    # --------------------------------------------------------
    # Dump matrix into a file
    # --------------------------------------------------------
    def save_lattice(self, filename):
        np.savetxt(filename, self.cluster_matrix(), fmt="%d")
