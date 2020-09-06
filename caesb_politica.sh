#!/bin/bash
#SBATCH -n 1 -N 1
#SBATCH --job-name=politica_96
#SBATCH --output=output/politica_96.out
#SBATCH --error=error/politica_96.err
#SBATCH --time=2:30:00
# #SBATCH --mail-user=andressa.giacomazzo@gmail.com
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=16
cd $SLURM_SUBMIT_DIR

time ./triangleSimulation\
	-T 16\
        -t 2086\
        -r 1000\
        -d /scratch/spec959/\
        -C -1\
	-m 0\
	-e 1\
	-s sample_solutions.csv\
        -O rof_tables_test_problem_intake_present/\
    	-U TestFiles/utilities_rdm.csv\
        -W TestFiles/water_sources_rdm.csv\
        -P TestFiles/policies_rdm.csv\
