#!/bin/bash
#SBATCH -n 16 -N 2
#SBATCH --cpus-per-task 2
#SBATCH --job-name=caesb_opt
#SBATCH --output=output/caesb_opt_1.out
#SBATCH --error=error/caesb_opt_1.err
#SBATCH --time=48:00:00
# #SBATCH --mail-user=andressa.giacomazzo@gmail.com
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=16
cd $SLURM_SUBMIT_DIR
time mpirun -np 16 ./triangleSimulation\
	-T 2\
	-t 2086\
	-r 1000\
	-d /scratch/spec959/\
	-C -1\
	-s sample_solutions.csv\
	-O rof_tables_test_problem_intake_present/\
	-e 0\
	-U TestFiles/utilities_rdm.csv\
	-W TestFiles/water_sources_rdm.csv\
	-P TestFiles/policies_rdm.csv\
	-o 50\
	-n 1000\
	-b true

