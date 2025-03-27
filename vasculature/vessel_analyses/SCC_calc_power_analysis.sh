#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify number of cores
#$ -pe omp 8
# Specify memory per core
#$ -l mem_per_core=8G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N power_analysis_reduced_difference

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -r "stats_power_analyses; exit"
