#!/bin/bash
#SBATCH -n 24 -N 8
#SBATCH --job-name=caesb_table
#SBATCH --output=output/borg_caesb_table.out
#SBATCH --error=error/borg_caesb_table.err
#SBATCH --time=12:00:00
#SBATCH --exclusive
# #SBATCH --mail-user=bct52@cornell.edu
# #SBATCH --mail-type=all
export OMP_NUM_THREADS=5
cd $SLURM_SUBMIT_DIR
module load valgrind/3.15.0
time mpirun -np 24 valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-file=valgrind_borg_caesb.out ./triangleSimulation\
	-T 5\
       -t 208\
       -r 10\
       -d /scratch/bct52/andressa/\
       -C -1\
       -O rof_tables_test_problem_intake_present/\
       -e 0\
       -U TestFiles/utilities_rdm.csv\
       -W TestFiles/water_sources_rdm.csv\
       -P TestFiles/policies_rdm.csv\
	-b true\
	-n 1000\
	-o 200
 
