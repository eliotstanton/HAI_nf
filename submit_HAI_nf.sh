#!/bin/bash
#SBATCH --time=8:00:01
#SBATCH --ntasks=1
#SBATCH --mem=2g
#SBATCH --tmp=2g

# --------------------------------------------------------------------------- #

# File:                 submit_HAI_nf.sh
# Date created:         03 June 2022
# Date last modified:   06 June 2022
# Author:               Eliot Stanton (eliot.stanton@state.mn.us)
# Description:          Slurm script for calling Nextflow QC pipeline

# --------------------------------------------------------------------------- #

module load nextflow/21.10.6

# Define arguments passed to the script from the command-line:
DIR_IN=$1

# Define filepaths used by this script:
BASE_DIR=/home/mdh/shared/software_modules/HAI_QC/0.nf
NEXTFLOW_CONFIG=$BASE_DIR/nextflow.config
DIR_OUT=HAI_nf_results

# Define help message:
HELP="submit_HAI_nf.sh [FASTQ_DIRECTORY]

    FASTQ_DIRECTORY: Directory containing FASTQ files

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
if [[ ! $1 ]]; then
        printf "ERROR: Incomplete arguments provided.\n\n"
        printf "$HELP\n\n"
        echo "Exiting."; exit

fi

# --------------------------------------------------------------------------- #

# Modify permissions for directory containing the FASTQ files:
chmod 770 -R $DIR_IN

# --------------------------------------------------------------------------- #

printf "HAI_QC.nf\n"
printf "\t-c $NEXTFLOW_CONFIG \ \n"
printf "\t--reads $DIR_IN \ \n"
printf "\t-resume\n"

HAI_QC.nf \
	-c $NEXTFLOW_CONFIG \
	--reads $DIR_IN \
	-resume

# --------------------------------------------------------------------------- #

# TODO: Iterate through DIR_OUT and change permissions:
chmod 770 -R $DIR_OUT
chmod 770 -R work
