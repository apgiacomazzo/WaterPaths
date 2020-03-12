#!/bin/bash
#SBATCH -n 1 -N 1 
#SBATCH --job-name=caesb_table_callgrind
#SBATCH --output=output/caesb_table_callgrind.out
#SBATCH --error=error/caesb_table_callgrind.err
#SBATCH --time=10:30:00
# #SBATCH --mail-user=bct52@cornell.edu
cd $SLURM_SUBMIT_DIR
module load valgrind/3.15.0
valgrind --tool=callgrind --callgrind-out-file=table_gen.callgrind ./triangleSimulation\
	-T 8\
       -t 586\
       -r 16\
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

valgrind --tool=callgrind --callgrind-out-file=use_table.callgrind ./triangleSimulation\
	-T 8\
       -t 586\
       -r 16\
       -d /scratch/bct52/andressa/\
       -C -1\
	-m 0\
	-s sample_solutions.csv\
       -O rof_tables_test_problem/\
       -e 0\
       -U TestFiles/utilities_rdm.csv\
       -W TestFiles/water_sources_rdm.csv\
       -P TestFiles/policies_rdm.csv
