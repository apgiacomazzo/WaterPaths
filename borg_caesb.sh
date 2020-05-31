#!/bin/bash
#SBATCH -n 48 -N 16
#SBATCH --job-name=caesb_table
#SBATCH --output=output/borg_caesb_table.out
#SBATCH --error=error/borg_caesb_table.err
#SBATCH --time=42:00:00
#SBATCH --exclusive
#SBATCH -c 5
# #SBATCH --mail-user=bct52@cornell.edu
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=5
cd $SLURM_SUBMIT_DIR

time mpirun -np 24 ./triangleSimulation\
	-T 5\
       -t 2086\
       -r 1000\
       -d /scratch/bct52/andressa/\
       -C -1\
       -O rof_tables_test_problem/\
       -e 1\
       -U TestFiles/utilities_rdm.csv\
       -W TestFiles/water_sources_rdm.csv\
       -P TestFiles/policies_rdm.csv\
	-b true\
	-n 75000\
	-o 1000
