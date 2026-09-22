# site_percolation_NZ.py

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
def _find(site, parent, dx, dy):
    """Find root and displacement of site relative to the root."""
    root = site
    x = 0
    y = 0

    while parent[root] >= 0:
        x += dx[root]
        y += dy[root]
        root = parent[root]

    # Path compression
    node = site
    px = 0
    py = 0

    while parent[node] >= 0:
        next_node = parent[node]

        old_dx = dx[node]
        old_dy = dy[node]

        parent[node] = root
        dx[node] = x - px
        dy[node] = y - py

        px += old_dx
        py += old_dy
        node = next_node

    return root, x, y


# ============================================================
# Numba-compileable add_site function. It is called outside Numba
# ============================================================

@njit(cache=True)
def _add_site(L, site, parent, dx, dy, occupied, wrap):
    """
    Occupy one site and update the weighted union-find structure.

    Returns
    -------
    new_size
        Size of the cluster containing the new site.
    erased_sizes
        Array containing the sizes of clusters that disappeared by merging.
    n_erased
        Number of entries in erased_sizes that are valid.
    wraps
        True if the resulting cluster wraps around the torus.
    """

    occupied[site] = True

    x = site % L
    y = site // L

    # Neighbours: site, dx, dy
    neighbours = np.empty((4, 3), dtype=np.int64)

    neighbours[0, 0] = (y * L + (x - 1) % L)
    neighbours[0, 1] = -1
    neighbours[0, 2] = 0

    neighbours[1, 0] = (y * L + (x + 1) % L)
    neighbours[1, 1] = 1
    neighbours[1, 2] = 0

    neighbours[2, 0] = (((y - 1) % L) * L + x)
    neighbours[2, 1] = 0
    neighbours[2, 2] = -1

    neighbours[3, 0] = (((y + 1) % L) * L + x)
    neighbours[3, 1] = 0
    neighbours[3, 2] = 1

    # The new site initially forms a cluster of size 1.
    parent[site] = -1
    dx[site] = 0
    dy[site] = 0



    erased_sizes = np.empty(4, dtype=np.int64)
    n_erased = 0

    root = site

    for k in range(4):

        neighbour = neighbours[k, 0]

        if not occupied[neighbour]:
            continue

        ex = neighbours[k, 1]
        ey = neighbours[k, 2]

        root_i, xi, yi = _find(site, parent, dx, dy)
        root_j, xj, yj = _find(neighbour, parent, dx, dy)

        if root_i == root_j:

            # Closed path around the torus.
            loop_x = xi + ex - xj
            loop_y = yi + ey - yj

            if loop_x != 0 or loop_y != 0:
                wrap[root_i] = True

            continue

        # Sizes before merging.
        size_i = -parent[root_i]
        size_j = -parent[root_j]


        erased_sizes[n_erased] = size_j
        n_erased += 1
        


        if size_i > size_j:
            # root_j becomes child of root_i
            parent[root_i] = -(size_i + size_j)

            # Displacement root_j -> root_i
            dx[root_j] = xi + ex - xj
            dy[root_j] = yi + ey - yj

            wrap[root_i] = wrap[root_i] or wrap[root_j]
            parent[root_j] = root_i

            root = root_i

        else:
            # root_i becomes child of root_j
            parent[root_j] = -(size_i + size_j)

            # Displacement root_i -> root_j
            dx[root_i] = xj - xi - ex
            dy[root_i] = yj - yi - ey

            wrap[root_j] = wrap[root_j] or wrap[root_i]
            parent[root_i] = root_j

            root = root_j

        new_size = size_i + size_j

    new_size = -parent[root]

    return new_size, erased_sizes, n_erased, wrap[root]


# ============================================================
# Public Python class (w API)
# ============================================================

class SitePercolation2D:

    def __init__(self, L, seed=None):

        self.L = L
        self.N = L * L

        # Random ordering of all sites.
        rng = np.random.default_rng(seed)
        self.permutation = rng.permutation(self.N)

        # Weighted union-find.
        # Negative value = cluster size.
        self.parent = np.full(self.N, -1, dtype=np.int64)

        # Displacement of each site relative to its parent.
        self.dx = np.zeros(self.N, dtype=np.int64)
        self.dy = np.zeros(self.N, dtype=np.int64)

        # Occupied sites.
        self.occupied = np.zeros(self.N, dtype=np.bool_)

        # Does the cluster rooted at a site wrap around the torus?
        self.wrap = np.zeros(self.N, dtype=np.bool_)


        # Number of occupied sites.
        self.n = 0

    @property
    def p(self):
        return self.n / self.N



    # --------------------------------------------------------
    # Add a site to matrix. Outputs new_size, erased_sizes, n_erased, wraps
    # --------------------------------------------------------

    def add_site(self):

        site = self.permutation[self.n]

        new_size, erased_sizes, n_erased, wraps = _add_site(self.L, site, self.parent, self.dx, self.dy, self.occupied, self.wrap)

        self.n += 1

        erased_sizes = tuple(erased_sizes[:n_erased])

        return new_size, erased_sizes, wraps


    # --------------------------------------------------------
    # Output matrix of labels
    # --------------------------------------------------------

    def cluster_matrix(self):

        labels = np.zeros(self.N, dtype=np.int64)

        for site in range(self.N):

            if self.occupied[site]:

                root, _, _ = _find(site, self.parent, self.dx, self.dy)

                labels[site] = -self.parent[root]

        return labels.reshape(self.L, self.L)


    # --------------------------------------------------------
    # Dump matrix into a file
    # --------------------------------------------------------
    def save_lattice(self, filename):
        np.savetxt(filename, self.cluster_matrix(), fmt="%d")


