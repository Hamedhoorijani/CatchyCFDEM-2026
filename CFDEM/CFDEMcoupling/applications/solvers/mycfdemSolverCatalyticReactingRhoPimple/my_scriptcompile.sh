#!/bin/bash -l
#
#PBS -N solver-compilation
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:59:59

source $HOME/setup_CFDEM-8_env-rome.sh

# Go to case directory
cd $PBS_O_WORKDIR

wclean
wmake
