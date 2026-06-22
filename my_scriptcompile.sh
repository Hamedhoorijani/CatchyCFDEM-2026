#!/bin/bash -l
#
#PBS -N rome-2023a-compilation
#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:59:59

source $HOME/setup_CFDEM-8_env-rome-foss-2023a.sh

# Go to case directory
cd $PBS_O_WORKDIR

./Allwclean
./Allwmake

cd CFDEM
source $VSC_DATA/catchyFOAM-rome-foss-2023a/CFDEM/CFDEMcoupling/etc/bashrc
bash $CFDEM_PROJECT_DIR/etc/compileCFDEMcoupling_all.sh
cfdemCompCFDEMall
cd CFDEM/CFDEMcoupling/src/lagrangian/cfdemParticleComp
wmake libso
