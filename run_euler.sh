#!/bin/bash

# Name of the excecutable
RUN=$1 

# Name of the template directory, which has 
# to have the executable and all the necessary input files. 
TEMP_DIR=$2

# Number of societies in the ensemble (Number of different simulations to be run in parallel)
NSC=$3

# The simulation time, which is to be asked of the cluster (-W param.), in format HH.MM 
TIME=$4


# Now we make a loop over 1..NSC and copy the template directory on each case and just run the simulation 
# in each of the new directories 

for ((i=1; i <= ${NSC}; i++))
do
	# Copy the template directory
	cp -r ./${TEMP_DIR} ./${TEMP_DIR}_run_$i
	
	# Go to the new directory and run the simulation
	cd ./${TEMP_DIR}_run_$i
	bsub -W ${TIME} -n 1 ./${RUN}
	#echo "-W ${TIME} -n 1 ./${RUN}"
	cd ../

done

