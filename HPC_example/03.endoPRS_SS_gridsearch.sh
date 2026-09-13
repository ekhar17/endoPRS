#!/bin/bash

#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=20GB
#SBATCH -t 3-
#SBATCH --array=1-25
#SBATCH --output=slurm_output/endoPRS_SS/%A_%a.txt

module load r

## $1 = thresh, $2 = thresh_endo, $3 = save_folder
## The array index selects which row of the 5x5 weight grid to fit.
Rscript 03.endoPRS_SS_gridsearch.R $1 $2 $SLURM_ARRAY_TASK_ID $3

## Submit every threshold pair, keeping the endophenotype threshold at least as stringent:
## sbatch 03.endoPRS_SS_gridsearch.sh 0.01 0.01 /path/to/results
## sbatch 03.endoPRS_SS_gridsearch.sh 0.01 1e-4 /path/to/results
## sbatch 03.endoPRS_SS_gridsearch.sh 0.01 1e-6 /path/to/results
## sbatch 03.endoPRS_SS_gridsearch.sh 1e-4 1e-4 /path/to/results
## sbatch 03.endoPRS_SS_gridsearch.sh 1e-4 1e-6 /path/to/results
## sbatch 03.endoPRS_SS_gridsearch.sh 1e-6 1e-6 /path/to/results

## Once all jobs are done, combine them:
## Rscript -e 'library(endoPRS); res <- select_endoPRS_ss("/path/to/results"); write.csv(res$beta, "endoPRS_ss_final_beta.csv", row.names = FALSE)'
