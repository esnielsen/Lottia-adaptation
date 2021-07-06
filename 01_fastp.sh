#!/bin/bash
#SBATCH --mail-user=esnielsen@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J fast-p
#SBATCH -p RM-shared
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=6
set -e
set -x
# To Run
# sbatch 01_fastp.sh

#Set up directory
cd /ocean/projects/deb200006p/enielsen/LGwork

# Call fastp package
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate fastp.filter

# Iterate over files in data folder
for file in $(ls -1 40_raw_data/*_R1*.gz)
do
    input_file=$(echo "$file" | perl -pe 's/_R1.*\.fastq.gz//')
    output_file=$(basename "$input_file")
    echo "Treating: $output_file"

    fastp -w 6 \
        -i "$input_file"_R1.fastq.gz \
        -I "$input_file"_R2.fastq.gz \
        -o 05_trimmed/"$output_file"_1.trimmed.fastq.gz \
        -O 05_trimmed/"$output_file"_2.trimmed.fastq.gz \
        -j 05_trimmed/01_reports/"$output_file".json \
        -h 05_trimmed/01_reports/"$output_file".html
done
