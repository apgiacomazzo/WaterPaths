#!/bin/bash
#SBATCH -n 1 -N 1 
#SBATCH --job-name=caesb_table
#SBATCH --output=output/caesb_table.out
#SBATCH --error=error/caesb_table.err
#SBATCH --time=02:30:00
# #SBATCH --mail-user=bct52@cornell.edu
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=16
cd $SLURM_SUBMIT_DIR
time ./triangleSimulation\
	-T 16\
       -t 2086\
       -r 1000\
       -d /scratch/bct52/andressa/\
       -C 1\
	-m 0\
	-s sample_solutions.csv\
       -O rof_tables_test_problem/\
       -e 0\
       -U TestFiles/utilities_rdm.csv\
       -W TestFiles/water_sources_rdm.csv\
       -P TestFiles/policies_rdm.csv\
	-p false

time ./triangleSimulation\
	-T 16\
       -t 2086\
       -r 1000\
       -d /scratch/bct52/andressa/\
       -C -1\
	-m 0\
	-s sample_solutions.csv\
       -O rof_tables_test_problem/\
       -e 0\
       -U TestFiles/utilities_rdm.csv\
       -W TestFiles/water_sources_rdm.csv\
       -P TestFiles/policies_rdm.csv
