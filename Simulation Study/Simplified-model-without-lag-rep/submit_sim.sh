#!/bin/sh

#SBATCH -J simulation
#SBATCH --array=1-20
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --mail-user=siyingma@uvic.ca
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm-%A.out
#SBATCH --account=def-lcowen

echo "The current directory is `pwd`"

module load r/4.4.0
export R_LIBS=/home/siyingma/projects/def-lcowen/siyingma/simplified

Rscript ./run_sim.R $SLURM_ARRAY_TASK_ID
