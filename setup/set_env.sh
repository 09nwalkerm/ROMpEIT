#!/bin/bash

echo "Setting up local environment..."
echo "Finding top of ROMEG tree..."

cd ..
export ROMEG=$PWD
echo " ----------------------------------------------------------------------"
echo " |                  Welcome to the ROMpEIT Repository                 |"
echo " |                                                                    |"
echo " | Please enter the path where you would like the data to be saved.   |"
echo " |                                                                    |"
echo " | Recommended options are either the top of the ROMpEIT repo tree in |"
echo " | your home directory, or another folder in which you would like the |"
echo " | results and intermediate results to be held in (e.g. a scratch     |"
echo " | storage pool).                                                     |"
echo " |                                                                    |"
echo " | Example: /path/to/scratch/storage/ROM_results                      |"
echo " |                                                                    |"
echo " | Alternatively leave blank to set the data path to the Results      |"
echo " | folder in the ROMpEIT repository.				    |"
echo " ----------------------------------------------------------------------"

read -e -p "Path (use tab for completion): " path

if [[ $path = "" ]]; then
	export ROMEG_DATA=$ROMEG/results
       	
	if [ -d ${ROMEG}/results ]; then
		echo "Using ROMpEIT/results folder for data storage."
	else
		echo "Making ROMpEIT/results folder for data storage."
		mkdir ${ROMEG}/results
	fi
else
	export ROMEG_DATA=$path
fi

if [ -d ${ROMEG_DATA}/logs ] || [ -d ${ROMEG_DATA}logs ]; then
	echo "Using /logs folder for logging."
	end=`echo $ROMEG_DATA | rev | cut -c 1`
	if [[ end = "/" ]]; then
		export ROMEG_LOGS=${ROMEG_DATA}logs
	else
		export ROMEG_LOGS=${ROMEG_DATA}/logs
	fi
else
	read -p "There is no logging folder in this directory, would you like to create one? (y/n): " ans
	if [[ $ans = "y" ]] || [[ $ans = "Y" ]] || [[ $ans = "yes" ]]; then
		end=`echo $ROMEG_DATA | rev | cut -c 1`
	        if [[ end = "/" ]]; then
			mkdir ${ROMEG_DATA}logs
			export ROMEG_LOGS=${ROMEG_DATA}logs
		else
			mkdir ${ROMEG_DATA}/logs
			export ROMEG_LOGS=${ROMEG_DATA}/logs
		fi
	fi
fi

