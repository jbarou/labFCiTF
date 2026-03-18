import numpy as np


#USE_NUMBA = False ### Compila el codi. Per fer proves, no cal. Quan feu les simualcions llargues, poseu a True.
USE_NUMBA = True ### Compila el codi. Per fer proves, no cal. Quan feu les simualcions llargues, poseu a True.



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


# ---------------------------------------------------------
# Numba compilable cluster growth kernel
# ---------------------------------------------------------

@njit(cache=True)
def grow_cluster_func(L, p, seed_x, seed_y):

    occupied = np.zeros((L, L), np.int8)

    # displacement arrays
    dx_arr = np.zeros((L, L), np.int32)
    dy_arr = np.zeros((L, L), np.int32)

    # queue (max possible cluster size L^2)
    qx = np.empty(L * L, np.int32)
    qy = np.empty(L * L, np.int32)

    head = 0
    tail = 0

    qx[tail] = seed_x
    qy[tail] = seed_y
    tail += 1

    occupied[seed_x, seed_y] = 1
    dx_arr[seed_x, seed_y] = 0
    dy_arr[seed_x, seed_y] = 0

    wraps = False

    # expands cluster using two pointers --(head*)-->(tail*)
    # propagates at head* and fills to tail* 
    while head < tail:

        x = qx[head]
        y = qy[head]
        head += 1

        dx0 = dx_arr[x, y]
        dy0 = dy_arr[x, y]

        # for the 4 neighbors
        for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:

            nx = (x + dx) % L
            ny = (y + dy) % L

            ddx = dx0 + dx
            ddy = dy0 + dy

            if occupied[nx, ny] == 1:

                if dx_arr[nx, ny] != ddx or dy_arr[nx, ny] != ddy:
                    wraps = True

            elif occupied[nx, ny] == 0:

                if np.random.random() >= p:
                    occupied[nx, ny] = -1
                    continue

                occupied[nx, ny] = 1
                dx_arr[nx, ny] = ddx
                dy_arr[nx, ny] = ddy
                occupied[nx,ny] = 1
                qx[tail] = nx
                qy[tail] = ny
                tail += 1


    return tail, wraps, occupied


# ---------------------------------------------------------
# Python wrapper class
# ---------------------------------------------------------

class SitePercolation2D:

    def __init__(self, L, p, seed=None):
        self.L = L
        self.p = p
        if seed is not None:
            np.random.seed(seed)

    # ---------------------------------------------------------
    # Exit function for 1 configuration. returns size, wraps, lattice
    # ---------------------------------------------------------

    def run_one(self, p=None):
        if p is not None:
           self.p = p
        size, wraps, lattice = grow_cluster_func(self.L, self.p, int(self.L/2), int(self.L/2))
        return size, wraps, lattice


    # ---------------------------------------------------------
    # Exit function for N configuration: statistics. Returns lists of sizes_fin, sizes_wra, n_wrapping
    # ---------------------------------------------------------

    def run_n(self, n_trials):

        sizes_fin = []
        sizes_wra = []
        n_wrapping = 0


        for k in range(n_trials):

            size, wraps, _ = grow_cluster_func(self.L, self.p,  int(self.L/2), int(self.L/2))

            if wraps:
                n_wrapping += 1
                sizes_wra.append(size)
            else:
                sizes_fin.append(size)

        return sizes_fin, sizes_wra, n_wrapping
