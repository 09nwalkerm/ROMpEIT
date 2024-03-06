#!/bin/bash

#SBATCH --job-name=CTRLleaf
#SBATCH --requeue
#SBACTH --begin=now
#SBATCH --partition=cubric-rocky8
#SBATCH -o ROMEG_%j.out
#SBATCH -e ROMEG_%j.err
#SBATCH --cpus-per-task=12

script=`echo $1 |awk '{split($1,a,"."); print a[1]; exit}'`

/opt/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -nosplash -r \
"addpath(genpath([getenv('ROMEG') '/Functions'])); \
addpath(genpath([getenv('ROMEG') '/scripts'])); \
addpath(genpath([getenv('ROMEG') '/external'])); \
$script; \
exit;"
