#!/bin/sh
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --time=3:00:00
#SBATCH --partition=batch
#SBATCH --licenses=common
#SBATCH --job-name=genomic_prediction_jyang_%A_%a
#SBATCH --mail-user=zhikaiyang911@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o /common/jyanglab/zhikaiyang/projects/GP_microbiome/slurm-log/genomic_prediction_jyang_%A_%a_out.txt
#SBATCH -e /common/jyanglab/zhikaiyang/projects/GP_microbiome/slurm-log/genomic_prediction_jyang_%A_%a_error.txt

#SBATCH --array=7-14

module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/GP_microbiome/largedata/genomic_prediction_jyang_data_seeds_vis.R $1 $2 $SLURM_ARRAY_TASK_ID
