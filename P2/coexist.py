

#### No es suficient, perque sempre hi queden mesures de "coexistencia" (possiblement occillant)
#### possibles solucions:
####    - fer passos de T més amples
####    - augmentar més el temps de termalitzacio
####    - complicar una mica el calcul, utilitzar thresholds o smoothcurves
####    - emprar variança en comptes de derivades en T


import glob
import numpy as np
import matplotlib.pyplot as plt


rows_out = []  # guardarem directament les columnes seleccionades

#files = sorted(glob.glob("thermo_NPT/avg_N125P*.dat"))
files = sorted(glob.glob("thermo_NPT_isoP/avg_N125P*.dat"))

for f in files:

    dataIn= np.loadtxt(f)
    col = dataIn[:,6]
    
    dcol = []

    for i in range(col.size-2):
      dcol.append(col[i+2]-col[i])

    # Índexs de la variació mínima i màxima
    i_min = np.argmin(dcol)
    i_max = np.argmax(dcol)

    # Les variacions corresponen al salt entre fila i i fila i+1
    row_min = dataIn[i_min]
    row_max = dataIn[i_max]

    diff = row_max - row_min
    ave  = 0.5*(row_max + row_min)

    # Escollim només les columnes 3,5,7 (índexs 2,4,6)
    out_row = [
        row_min[2], row_min[4], row_min[6],
        row_max[2], row_max[4], row_max[6],
        ave[2], ave[4], ave[6],
        diff[2], diff[4], diff[6]
    ]

    rows_out.append(out_row)


# Convertim a array 2D un cop al final
dataOut = np.array(rows_out)

np.savetxt("coex.out", dataOut, fmt="%.6f")
