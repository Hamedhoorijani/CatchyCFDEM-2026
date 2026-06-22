#!/bin/bash -l
#
#PBS -N liggghtsParallelTest
#PBS -l nodes=1:ppn=4
#PBS -l walltime=71:59:59

source $HOME/setup_CFDEM_8_shinx.sh

cd $PBS_O_WORKDIR


# run liggghts init
mympirun liggghts -in in.heatGran

