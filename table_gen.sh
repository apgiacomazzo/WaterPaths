#!/bin/bash
#SBATCH -n 1 -N 1 
#SBATCH --job-name=caesb_table
#SBATCH --output=output/caesb_table.out
#SBATCH --error=error/caesb_table.err
#SBATCH --time=5:30:00
# #SBATCH --mail-user=apgiacomazzo@gmail.com
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=16
cd $SLURM_SUBMIT_DIR
time ./triangleSimulation\
	-T 16\
        -t 2086\
        -r 1000\
        -d /scratch/spec959/\
        -C 1\
	-m 0\
	-s sample_solutions.csv\
        -O rof_tables_test_problem_intake_present/\
        -e 0\
        -U TestFiles/utilities_rdm.csv\
        -W TestFiles/water_sources_rdm.csv\
        -P TestFiles/policies_rdm.csv\
	-p false

time ./triangleSimulation\
	-T 16\
        -t 2086\
        -r 1000\
        -d /scratch/spec959/\
        -C -1\
	-m 0\
	-s sample_solutions.csv\
        -O rof_tables_test_problem_intake_present/\
        -e 0\
        -U TestFiles/utilities_rdm.csv\
        -W TestFiles/water_sources_rdm.csv\
        -P TestFiles/policies_rdm.csv
