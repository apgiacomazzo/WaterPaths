#!/bin/bash
#SBATCH -n 1 -N 1
#SBATCH --job-name=calibrate
#SBATCH --output=output/calibrate.out
#SBATCH --error=error/calibrate.err
#SBATCH --time=00:16:00
cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=16
N_REALIZATIONS=64
SOLUTION=0
DATA_DIR="/home/fs02/pmr82_0001/bct52/andressa/WaterPaths/"
./triangleSimulation -T ${OMP_NUM_THREADS}\
       	-t 1045\
       	-r ${N_REALIZATIONS}\
       	-d ${DATA_DIR}\
       	-C 0\
       	-s sample_solutions.csv\
       	-e 0\
        -m ${SOLUTION}\
        -p false	
# ./triangleSimulation -T ${OMP_NUM_THREADS}\
#        	-t 2344\
#        	-r ${N_REALIZATIONS}\
#        	-d ${DATA_DIR}\
#        	-C -1\
# 	-O rof_tables_test_problem\
#        	-s reference_final.reference\
#        	-e 0\
# 	-f 0\
# 	-l 209
# 	-U TestFiles/rdm_utilities_test_problem_opt.csv\
# 	-W TestFiles/rdm_water_sources_test_problem_opt.csv\
#         -P TestFiles/rdm_dmp_test_problem_opt.csv\
#         -p false	
