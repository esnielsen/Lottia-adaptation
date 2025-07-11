                                #!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J pbscan
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=40
set -e
set -x
# To Run
# sbatch pbscan.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/Outliers/PBScan/ca_tests


./pbscan -likes ALL.likes -pop1 CA.n.pbscan.txt -pop2 CA.c.pbscan.txt -pop3 MEX$

# the program takes the following input files
# - likes = the beagle file from angsd (need to change suffix to .likes)
# - pops 1 & 2 = list of individuals in each of the 2 populations you want to c$
# - pop3 is the list of individuals in the outgroup
# - out =output name
# - ms= neutral data generated with msms program
# - min = minimum number of individuals need to compare allele freqs
# - maf = min minor allele freq

