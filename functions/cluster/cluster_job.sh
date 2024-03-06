#!/bin/bash

#SBATCH --no-requeue
#SBACTH --begin=now
#SBATCH --partition=cubric-rocky8
#SBATCH --exclude=c4b[8,17],c5b[2-6,8-10,12-17],c6b6,c6b7
###SBATCH --cpus-per-task 2
###SBATCH --ntasks-per-node=10


cd $ROMEG_TOP

echo This is running on node `hostname`

/opt/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -nosplash -r \
"addpath(genpath([getenv('ROMEG') '/Functions'])); \
addpath(genpath([getenv('ROMEG') '/external'])); \
cluster_job_$1(); \
exit;"
