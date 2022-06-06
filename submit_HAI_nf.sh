#!/bin/bash
#SBATCH --time=8:00:01
#SBATCH --ntasks=1
#SBATCH --mem=2g
#SBATCH --tmp=2g

# --------------------------------------------------------------------------- #

# File:                 submit_HAI.sh
# Date created:         03 June 2022
# Date last modified:   06 June 2022
# Author:               Eliot Stanton (eliot.stanton@state.mn.us)
# Description:          Slurm script for calling Nextflow QC pipeline

# --------------------------------------------------------------------------- #

module load nextflow

# Define arguments passed to the script from the command-line:
DIR_IN=$1
DIR_OUT=$2

# Define variable for error reporting:
VAR_ERROR=0

# Define help message:
HELP="submit_HAI_nf.sh [FASTQ_DIRECTORY] [OUTPUT_DIRECTORY]

    FASTQ_DIRECTORY: Directory containing FASTQ files
    ACCESSION: Accession number of sample to be analyzed (YYYYEL-####)
    OUTPUT_DIRECTORY: A subdirectory containing results will be created here

    FASTQ files should be formatted with the following formatting:
	 YYYYEL-####_R1.fastq.gz and YYYY-####_R2.fastq.gz.

    Files in OUTPUT_DIRECTORY/ACCESSION include:
	ACCESSION_filtered_R1.fastq.gz
	ACCESSION_filtered_R2.fastq.gz
	ACCESSION_trimmed.fa
	ACCESSION_kraken2_report.txt
	ACCESSION_stats.txt"

# --------------------------------------------------------------------------- #

# Handle not enough arguments being supplied to the script:
if [[ ! $2 ]]; then
        printf "ERROR: Incomplete arguments provided.\n\n"
        printf "$HELP\n\n"
        echo "Exiting."; exit

fi

# --------------------------------------------------------------------------- #

# Modify permissions for directory to contain processed files:
chmod 770 -R $DIR_OUT

# Modify permissions for directory containing the FASTQ files:
chmod 770 -R $DIR_IN

# --------------------------------------------------------------------------- #

printf "HAI_QC.nf\n"
printf "\t-c nextflow.config \ \n"
printf "\t--reads $DIR_IN \ \n"
printf "\t-resume\n"

./HAI_QC.nf \
	-c nextflow.config \
	--reads $DIR_IN \
	-resume

# --------------------------------------------------------------------------- #

# Change permissions recursively for all the new files within $DIR_OUT2:
chmod -R 770 $DIR_OUT
