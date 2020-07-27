#!/bin/bash
#SBATCH -n 24 -N 8
#SBATCH --job-name=caesb_semFC_1
#SBATCH --output=output/caesb_semFC_1.out
#SBATCH --error=error/caesb_semFC_1.err
#SBATCH --time=60:00:00
#SBATCH --exclusive
#SBATCH -c 5
# #SBATCH --mail-user=andressa.giacomazzo@gmail.com
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=5
cd $SLURM_SUBMIT_DIR

time mpirun -np 24 ./triangleSimulation\
	-T 5\
        -t 2086\
        -r 1000\
        -d /scratch/spec959/\
        -C -1\
        -O rof_tables_test_problem_intake_present/\
        -e 1\
        -U TestFiles/utilities_rdm.csv\
        -W TestFiles/water_sources_rdm.csv\
        -P TestFiles/policies_rdm.csv\
	-b true\
	-n 35000\
	-o 1000
