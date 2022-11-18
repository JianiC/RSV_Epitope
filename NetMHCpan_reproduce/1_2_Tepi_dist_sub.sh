#!/bin/bash

#SBATCH --job-name=B0702
#SBATCH --partition=bahl_salv_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --cpus-per-task=4
#SBATCH --time=150:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

cd $SLURM_SUBMIT_DIR

module load Miniconda2/4.7.10
source activate /home/jc54391/Tepi_env

python Tepitope_dis_allele.py -i netMHCpan_out/HLAB4403.xls -o netMHCpan_out/HLAB4403_seqdist.csv
conda deactivate
