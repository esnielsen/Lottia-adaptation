#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J LDannot
#SBATCH -p RM-shared
#SBATCH -t 8:00:00
#SBATCH --ntasks-per-node=1
set -e
set -x
# To Run
# sbatch LDannot.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork/LDannot/LD-annot

# Call package
module load anaconda3

#set path to plink
export PATH=/jet/home/enielsen:$PATH

#run LDannot
python3 LD-annot0.4.py ALL_2_no_chr.vcf Lottia_gigantea.55.edit.gff harv.outs.LDannot.txt gene 0.9 harv.out.annot

# from https://github.com/ArnaudDroitLab/LD-annot: 
# python3 LD-annot0.4.py PathToSnpFiles ex_annotations.gff ex_candidate_SNPs.txt type thr output SNP_Map

# "type" is the feature (mRNA, CDS, gene)
# "thr" is r2 threshold 
