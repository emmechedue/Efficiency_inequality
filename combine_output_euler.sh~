#!/bin/bash

# Name of the template dir (the base name for the simulation) 
TEMP_DIR=$1



# Make the output dir 
DIR=./${TEMP_DIR}_combined
mkdir ./${DIR}


# Make the output files and headers 
cd ./${TEMP_DIR}_run_1/
for file in *.txt 
do 
	sed -n 1,4p ${file} > ../${DIR}/${file}

done
cd ../


# Combine the outputs
for DSC in ${TEMP_DIR}_run_*
do

	# Get info for each file form society
	cd ./${DSC}/
	for file in *.txt 
	do 
		sed -n 5p ${file} >> ../${DIR}/${file}

	done
	cd ../

done
