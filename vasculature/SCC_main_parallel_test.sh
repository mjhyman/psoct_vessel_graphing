#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=1:00:00

# Name of job
#$ -N test_array

# Combine output/error files into single file
#$ -j y

### Declare array job (create new job for each)
# There are 17 subjects and 3 gaussian sigma arrays
# (small, medium, large) for a total of 51 jobs. 
#$ -t 1-9

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load matlab/2022b

### Attempt simple SGE_TASK_ID syntax
matlab -nodisplay -batch scc_parallel_test $SGE_TASK_ID

### Attempts at iterating a list of two variables
# Retrieve contents of line number $SGE_TASK_ID from the text file subject_sigma_index_list
# params=`sed -n "${SGE_TASK_ID} p" subject_sigma_index_list.txt`
# paramsArray=($params)

# Assign first element of line to subject ID index
# subid_idx=${paramsArray[0]}

# Assign second element of line to gaussian sigma array index
# gauss_idx=${paramsArray[1]}

# The next line did not pass in the arguments subid_idx and gauss_idx
# matlab -nodisplay -batch "scc_parallel_test" $subid_idx $gauss_idx

# Attempt a different syntax
# matlab -nodisplay -batch "subid_idx='$subid_idx'; gauss_idx='$gauss_idx'; scc_parallel_test"

