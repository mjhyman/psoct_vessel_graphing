#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 28 cores and at least 384 GB of RAM.
# Specify number of cores
#$ -pe omp 28
# Specify memory per core
#$ -l mem_per_core=13G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=24:00:00

# Name of job
#$ -N gen_graph

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -singleCompThread -r "gen_graph_main; exit"
