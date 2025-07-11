#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J baypass.env
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=5
set -e
set -x
# To Run
# sbatch baypass.env.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/Outliers/03_baypass

NB_CPU=5


baypass_2.31/sources/g_baypass -npop 19 \
 -gfile by_pop_0.05_pctind0.5_maxdepth15.mafs.pruned.baypass \
 -omegafile ALL_2_r2_mat_omega.out \
 -efile bpass.env.NEW.txt  -scalecov \
 -outprefix ALL_2_r2.NEW.new.env.output -nthreads $NB_CPU
