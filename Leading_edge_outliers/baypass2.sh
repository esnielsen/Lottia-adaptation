#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J baypass2
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=5
set -e
set -x
# To Run
# sbatch baypass2.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/Outliers/03_baypass


baypass_2.31/sources/g_baypass -npop 2 -gfile G.CA.n.c_sim_pod -outprefix CA.n.c_r2.sim -npilot 100 -nthreads 5
