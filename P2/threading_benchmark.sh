#!/bin/bash

LAMMPS=lmp
#INPUT=in.threading.lj
#INPUT=in_v2.lj
INPUT=OpenMP_test.lj

export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "MPI OMP TotalCores Time(s)"

for MPI in 1 2 4; do
  for OMP in 1 2 4; do

    TOTAL=$((MPI * OMP))
    if [ $TOTAL -gt 8 ]; then
      continue
    fi

    export OMP_NUM_THREADS=$OMP

    echo "Running MPI=$MPI OMP=$OMP"

    mpirun -np $MPI $LAMMPS -sf omp -pk omp $OMP -in $INPUT \
      | tee out_mpi${MPI}_omp${OMP}.log

    # Extract total wall time
    grep "Loop time of" out_mpi${MPI}_omp${OMP}.log
    SECS=$(grep "Loop time of" out_mpi${MPI}_omp${OMP}.log | awk '{print $4}')

    echo "$MPI $OMP $TOTAL | $SECS" >> scaling_results.dat

  done
done

