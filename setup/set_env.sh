#!/bin/bash

echo "Setting up local environment..."
echo "Finding top of ROMEG tree..."

cd ..
export ROMEG=$PWD
echo " ----------------------------------------------------------------------"
echo " |                  Welcome to the ROMEG Repository                   |"
echo " |                                                                    |"
echo " | Please enter the path where you would like the data to be saved.   |"
echo " |                                                                    |"
echo " | Recommended options are either the top of the ROMEG repo tree in   |"
echo " | your home directory, or another folder in which you would like the |"
echo " | results and intermediate results to be held in (e.g. a scratch     |"
echo " | storage pool).                                                     |"
echo " |                                                                    |"
echo " | Example: /cubric/data/c1616132/ROMEG_data/                         |"
echo " |                                                                    |"
echo " | Alternatively leave blank to set the data path to the Results      |"
echo " | folder in the ROMEG repository.                                    |"
echo " ----------------------------------------------------------------------"

read -e -p "Path (use tab for completion): " path

if [[ $path = "" ]]; then
	export ROMEG_DATA=$ROMEG/Results 
else
	export ROMEG_DATA=$path
fi

if [ -d ${path}/logs ] || [ -d ${path}logs ]; then
	echo "Using /logs folder for logging."
else
	read -p "There is no logging folder in this directory, would you like to create one? (y/n): " ans
	if [[ $ans = "y" ]] || [[ $ans = "Y" ]] || [[ $ans = "yes" ]]; then
		end=`echo $path | rev | cut -c 1`
	        if [[ end = "/" ]]; then
			mkdir ${path}logs
		else
			mkdir ${path}/logs
		fi
	fi
fi

