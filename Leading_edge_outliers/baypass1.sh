#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J baypass1
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=5
set -e
set -x
# To Run
# sbatch baypass1.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/Outliers/03_baypass


#baypass_2.31/sources/g_baypass -npop 2 -gfile p.CA.n.c.baypass -outprefix ALL_2_r2 -npilot 100 -nthreads 5

#to get C2 constrast values 
baypass_2.31/sources/g_baypass -npop 15 -gfile by_pop_0.05_pctind0.5_maxdepth15.mafs.ALL_2.L.baypass -contrastfile harv_press_bpass.contrast.txt -outprefix harv_press_L_C2 -npilot 100 -nthreads 5

