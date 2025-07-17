#!/bin/bash

#SBATCH --job-name=Apr15_ddd-dmu--4.5_MUgb--6.0_gdrug--1_1
#SBATCH --output=Apr15_ddd-dmu--4.5_MUgb--6.0_gdrug--1_1.txt

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute,hagan-compute-short,hagan-compute-long
##SBATCH --nodelist=compute-9-3
#SBATCH --qos=medium

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --chdir=SM_high-DMU--4.5_MUGB--6.0_GDRUG--1

#module load GSL/2.5
module load share_modules/GSL/2.5 
echo "Start 10486856 @ `date`"
scp -r ../source seed-1 ;
sleep 10
cd seed-1 ;
cd source ;
make ;
sleep 10 
cd .. ;
# 0.404 0.247 -0.146 0 -0.663 -.8
./source/assemble 10486856 4200.000 40.000 800.000 0.240 0.480 -9.81 -11.5 0.002 -4.5 0.100 -1 -15 0.0000002 0.3 0.1 -0.1 0 -0.6 -.95
sleep 1
rm -rf source
echo "End 10486856 @ `date`"
