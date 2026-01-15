#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --time=336:00:00
#SBATCH --mail-user=kvo020@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu
module load anaconda
conda create --Khoi_env python=3.11
conda activate Khoi_env


# Capturing the first passed argument
#simulation_id=$1

# Running your simulation program with the captured argument
python3 driver_fungalGrowth_singleNutrient.py

conda deactivate