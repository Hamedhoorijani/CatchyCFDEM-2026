#!/bin/bash -l

#SBATCH --job-name=2D-ETFB-60v				#Job name
#SBATCH --partition=HPC_PARTITION_NAME			#Partition
#SBATCH --account=HPC_PROJECT_NAME	 		#Project Name
#SBATCH --nodes=1					#Number of nodes
#SBATCH --ntasks=16					#Number of cores
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00				#Time:day-hours:minutes:seconds
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set +e

# Use the submission directory unless a case directory is supplied.
CASE_DIR="${1:-${SLURM_SUBMIT_DIR}}"

source $HOME/catchyCFDEM.sh

cd "$CASE_DIR"

cd DEM 
srun liggghts -in in.liggghts_init

cd ../CFD
decomposePar -force > log.decomposePar
srun ETFBCfdemSolverRhoPimple -parallel 
