#!/bin/bash

#SBATCH --job-name=hpp
#SBATCH --time=17:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250m
#SBATCH	--array=1-900
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/hpp/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/hpp/Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/R/02_run_sims_hpp.R /nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/hpp/Rout/phat_$SLURM_ARRAY_TASK_ID.Rout