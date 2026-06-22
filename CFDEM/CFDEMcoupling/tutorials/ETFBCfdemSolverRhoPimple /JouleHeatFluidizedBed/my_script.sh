#!/bin/bash -l
#
#PBS -N ETFB-CG2
#PBS -l nodes=1:ppn=16
#PBS -l walltime=71:59:59
#PBS -A 2025_042

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
#echo "BEFORE SOLVER"
#mpirun -np $PBS_NP ETFBCfdemSolverRhoPimple -parallel > log.ETFBCfdemSolverRhoPimple
mympirun --output=log.ETFBCfdemSolverRhoPimple  ETFBCfdemSolverRhoPimple -parallel
