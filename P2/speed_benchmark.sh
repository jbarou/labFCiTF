#!/bin/bash

LAMMPS=lmp
INPUT=inSpeedTest.lj

export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "MPI OMP TotalCores Time(s)"


mkdir testSpeed

rm  testSpeed/threading_speeds.dat

echo "#MPI OMP TOTAL SPEED(timesteps/s) " >> testSpeed/threading_speeds.dat

for MPI in 1 2 4 8; do
  for OMP in 1 2 4 8; do

    TOTAL=$((MPI * OMP))
    if [ $TOTAL -le 8 ]; then

        export OMP_NUM_THREADS=$OMP
        echo "================================================================================"
        echo "Running MPI=$MPI OMP=$OMP"
        echo "================================================================================"
  
        mpirun -np $MPI \
             --bind-to core \
             --map-by core \
             $LAMMPS -sf omp -pk omp $OMP -in $INPUT \
          | tee testSpeed/out_mpi${MPI}_omp${OMP}.log

        sleep 0.5
        # Extract performace
        grep "Performance:" testSpeed/out_mpi${MPI}_omp${OMP}.log
        SPEED=$(grep "Performance:" testSpeed/out_mpi${MPI}_omp${OMP}.log | awk '{print $4}')

        echo "$MPI     $OMP    $TOTAL      $SPEED" >> testSpeed/threading_speeds.dat

    fi

  done
  sleep .1
done

