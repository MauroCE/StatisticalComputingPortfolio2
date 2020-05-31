#!/bin/bash
#
#
#PBS -l nodes=1:ppn=10,walltime=1:00:00



# Define working directory
export WORK_DIR=$HOME/cpp_job

# Define executable
export EXE1="g++ -fopenmp simplecpp.cpp -o simplecpp"
export OMP_NUM_THREADS=10
export EXE3="./simplecpp"

# Add R module
module add languages/gcc-9.1.0

# Change into working directory
cd $WORK_DIR

# Do some stuff
echo JOB ID: $PBS_JOBID
echo Working Directory: `pwd`
echo Start Time: `date`

# Execute code
$EXE1
$EXE2
$EXE3 

echo End Time: `date`



sleep 20