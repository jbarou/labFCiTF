import numpy as np
from collections import deque


class SitePercolation2D:
    """
    Leath (epidemic) site percolation on a 2D torus.
    """

    def __init__(self, L, p, seed=None):
        self.L = L
        self.p = p
        self.seed = seed or np.random.default_rng()
        np.random.seed(seed)

    # ---------------------------------------------------------
    def grow_cluster(self, seed):
        """
        Grow a single cluster starting from `seed`.

        Returns
        -------
        size : int
            Cluster size
        wraps : bool
            True if the cluster wraps the torus
        dist : dict[(int,int) -> (int,int)]
            Map of occupied sites to displacement vectors
        """

        occupied = np.zeros((self.L, self.L), dtype=int)

        L = self.L
        p = self.p

        queue = deque([seed])

        occupied[seed]=1


        # dist[(x,y)] = (dx, dy) displacement from seed
        dist  = {seed: (0, 0)}
        wraps = False

        while queue:
            x, y = queue.popleft()
            dx0, dy0 = dist[(x, y)]

            for dx, dy in [(1,0), (-1,0), (0,1), (0,-1)]:

                nx = (x + dx) % L
                ny = (y + dy) % L

                ddx = dx0 + dx
                ddy = dy0 + dy

                if occupied[nx, ny] == 1:
                    old_dx, old_dy = dist[(nx, ny)]
                    if old_dx != ddx or old_dy != ddy:
                        wraps = True
                # Bernoulli trial
                elif occupied[nx, ny] == 0:
                    if np.random.rand() >= p: # not connected (but visited!!)
                        occupied[nx, ny] = -1
                        continue
                    occupied[nx, ny] = 1
                    dist[(nx, ny)] = (ddx, ddy)
                    queue.append((nx, ny))

        #size = occupied.count(1)
        size =  np.sum(occupied == 1)
        return size,wraps, occupied

    # ---------------------------------------------------------
    def run_one(self, seed):
        """
        Grow a cluster and return its configuration matrix.

        Parameters
        ----------
        seed : tuple[int, int]

        Returns
        -------
        C : ndarray (L, L)
            Binary matrix of the cluster
        size : int
            Cluster size
        wraps : bool
            Whether the cluster wraps
        """

        size, wraps, lattice = self.grow_cluster(seed)


        return size, wraps, lattice

    # ---------------------------------------------------------
    def run_n(self, n_trials):
        """
        Run multiple independent cluster growths.

        Returns
        -------
        sizes : list[int]
            Cluster sizes
        n_wrapping : int
            Number of wrapping clusters
        """

        sizes_fin = []
        sizes_wra = []
        n_wrapping = 0

        for _ in range(n_trials):
            seed = (int(self.L/2),int(self.L/2))
            size, wraps, _ = self.grow_cluster(seed)

            if wraps:
                n_wrapping += 1
                sizes_wra.append(size)
            else: 
                sizes_fin.append(size)
            
            

        return sizes_fin, sizes_wra, n_wrapping


