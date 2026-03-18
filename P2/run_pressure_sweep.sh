#!/bin/bash

LAMMPS=lmp
MPI=2
OMP=2

export OMP_NUM_THREADS=$OMP

INPUT=in_NPT_isobars.lj

PRESS=(0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15)

for P in "${PRESS[@]}"; do
    echo "Running P = $P"
     mpirun -np $MPI $LAMMPS -var P $P -sf omp -pk omp $OMP -in $INPUT 
done
