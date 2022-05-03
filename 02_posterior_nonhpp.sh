#!/bin/bash

#SBATCH --job-name=sample_nonhpp
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500m
#SBATCH	--array=[262,325]
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/non_hpp/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/non_hpp/Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/R/02_run_sims_non_hpp.R /nas/longleaf/home/ethanalt/Projects/Paper3/extensive_sims/Results/sims_raw/non_hpp/Rout/phat_$SLURM_ARRAY_TASK_ID.Rout