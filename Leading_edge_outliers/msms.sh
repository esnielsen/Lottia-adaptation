#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J msms
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=1
set -e
set -x
# To Run
# sbatch msms.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/Outliers/PBScan

# make sure msms3.2rc-b163.jar is in your working dir

java -jar msms3.2rc-b163.jar -ms  449 1000 -t 0.64