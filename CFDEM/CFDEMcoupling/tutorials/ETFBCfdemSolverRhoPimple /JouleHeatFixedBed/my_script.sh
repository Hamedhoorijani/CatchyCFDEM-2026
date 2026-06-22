#!/bin/bash -l
#
#PBS -N efixedbed
#PBS -l nodes=1:ppn=2
#PBS -l walltime=71:59:59
#PBS -l mem=32gb


source $HOME/setup_CFDEM-8-dev_env-rome-foss-2023a.sh

# -------------------------
# Go to case directory
cd $PBS_O_WORKDIR

# run liggghts init
#cd DEM
#mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init
#cd ..

cd CFD
#cleanCase
#blockMesh > log.blockMesh
decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemSolverRhoPimple -parallel > log.cfdemSolverRhoPimple
#mympirun --output=log.cfdemSolverRhoPimple  cfdemSolverRhoPimple -parallel

