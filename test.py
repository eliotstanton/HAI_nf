#!/usr/bin/python3

import sys
import subprocess
#import requests
import pandas as pd
import os
import re

# --------------------------------------------------------------------------- #

# File name:		test.py
# Date created:		31 May 2022
# Date modified:	31 MAY 2022
# Author:		Eliot Stanton (eliot.stanton@state.mn.us)
# Description:		This script pulls genome deliverables from sequencing
#			QC pipeline for NCBI metadata.

# --------------------------------------------------------------------------- #

# Variables imported by this script:
var_accession	= sys.argv[1]
dir_in		= sys.argv[2]	
dir_out		= sys.argv[3]

#print (var_accession, dir_in)

# Filepaths used by this script:
file_kraken2	= "{}/{}/{}_kraken2_report.txt".format(dir_out, var_accession, var_accession)
file_FASTQ1	= "{}/{}/{}_filtered_R1.fastq".format(dir_out, var_accession, var_accession)
file_FASTQ2	= "{}/{}_R1_001.fastq.gz".format(dir_in, var_accession)
file_quast	= "{}/{}/report.txt".format(dir_out, var_accession )
#print(file_FASTQ1, "\n", file_FASTQ2)

# Data structures used by this script:


# Variables used by this script:
var_ori_reads		= 1
var_html		= "https://www.ncbi.nlm.nih.gov/genome/?term="
var_avg_length		= 0
var_length		= 0
var_GC			= 0
var_num_reads		= 0
var_coverage		= 0
var_dropped		= 0
var_genus		= "unknown"
var_species		= "unknown"
var_align_genus		= 0
var_align_species	= 0
var_ratio		= 0
var_date		= "date +%F"
var_date		= subprocess.getoutput(var_date)

# --------------------------------------------------------------------------- #

# Determine the number of original paired reads:
var_ori_reads	= "echo $(zcat {} | wc -l) /4 | bc".format(file_FASTQ2)
var_ori_reads	= subprocess.getoutput(var_ori_reads)
#print (var_ori_reads)

# --------------------------------------------------------------------------- #

# Define data list:
data_list = []
data_genus = []
data_species = []

# Read through file_kraken2:
with open(file_kraken2) as file_read:
	for line in file_read:
		# Remove leading and trailing spaces:
		line = line.strip()

		# Split string into a list using tabs:
		split_line = line.split('\t')

		# Get unclassified reads result (denoted by 'unclassified') and append to data_list:
		if split_line[5] == 'unclassified':
			data_list.append(split_line)

		# Get species results (denoted by 'S') and append to data
		if split_line[3] == 'S':
			data_list.append(split_line)
			data_species.append(split_line)

		# Get species results (denoted by 'G') and append to data
		if split_line[3] == 'G':
			data_list.append(split_line)
			data_genus.append(split_line)

# Convert data lists to data frames:
data_list_df    = pd.DataFrame(data_list, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
data_genus_df   = pd.DataFrame(data_genus, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
data_species_df = pd.DataFrame(data_species, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])

# Remove left leading spaces from the Name column
data_list_df['Name']    = data_list_df['Name'].str.lstrip()
data_genus_df['Name']   = data_genus_df['Name'].str.lstrip()
data_species_df['Name'] = data_species_df['Name'].str.lstrip()

# Sort dataframes by percentages:
data_list_df    = data_list_df.sort_values(by=['Percentage'], ascending=False)
data_genus_df   = data_genus_df.sort_values(by=['Percentage'], ascending=False)
data_species_df = data_species_df.sort_values(by=['Percentage'], ascending=False)

		# Confine to Name and Percentage columns:
list_df         = data_list_df[['Name', 'Percentage']]
genus_df        = data_genus_df[['Name', 'Percentage']]
species_df      = data_species_df[['Name', 'Percentage']]

# Get first match from each data frame:
list_df         = list_df.head(1)
genus_df        = genus_df.head(1)
species_df      = species_df.head(1)

# Pull values:
list_name	= list_df['Name'].iloc[0]
var_genus	= genus_df['Name'].iloc[0]
var_species	= species_df['Name'].iloc[0]

list_perc       	= list_df['Percentage'].iloc[0]
var_align_genus		= genus_df['Percentage'].iloc[0]
var_align_species	= species_df['Percentage'].iloc[0]

var_species    = var_species.split(' ')

var_species    = var_species[-1]

#print (f"{var_genus} {var_align_genus}% {var_species} {var_align_species}%")

# --------------------------------------------------------------------------- #

# Scrape average species genome length from NCBI:
var_html = var_html + var_genus + "+" + var_species
file_html	= var_accession + ".html"
wget	 = "wget --output-document {} -q {}.html".format(var_accession,var_html)
os.system(wget)

with open(file_html) as file_read:
	for line in file_read:
		search	= re.search("median total length",line)
		if search:
			sline	= re.split('\>|\<',line)
			for i in sline:
				search  = re.search("median total length",i)
				if search:
					j = re.split(': ', i) 
					var_avg_length	= j[-1]
					var_avg_length	= float(var_avg_length)
					var_avg_length	*= 1000000
					var_avg_length	= int(var_avg_length)

# --------------------------------------------------------------------------- #

# Iterate through Quast results:
with open(file_quast) as file_read:
	for line in file_read:

		# Total assembly length:
		search_total_length	= re.search("^Total length \(>= 0 bp\)",line)
		if search_total_length:
			var_length	= re.split(' ', line)
			var_length	= var_length[10]
#			print (var_length)

		# Number of contigs:
		search_contigs		= re.search("# contigs  ", line)
		if search_contigs:
			var_contigs	= re.split (' ', line)
			var_contigs	= var_contigs[20]
#			print (var_contigs)

		# GC content:
		search_GC		= re.search("GC ", line)
		if search_GC:
			var_GC		= re.split (' ', line)
			var_GC		= var_GC[23]
#			print (var_GC)

		# Post-processed read number:
		search_left		= re.search("# left", line)
		if search_left:
			var_num_reads	= re.split (' ', line)
			var_num_reads	= var_num_reads[23]
#			print (var_num_reads)

		# Read depth:
		search_coverage		= re.search("Avg. coverage depth", line)
		if search_coverage:
			var_coverage	= re.split(' ', line)
			var_coverage	= var_coverage[11]
#			print (var_coverage)

# --------------------------------------------------------------------------- #

# Valculate dropped read percentage and round to two decimals:
var_num_reads	= int(var_num_reads)
var_ori_reads	= int(var_ori_reads)
var_dropped	= 1-(var_num_reads/var_ori_reads)
var_dropped	= round(var_dropped, 2)

# Determine genome ratio and round to two decimals:
var_length	= int(var_length)
var_avg_length	= int(var_avg_length)
if ( var_avg_length >0 ):
	var_ratio	= var_length/var_avg_length
var_ratio	= round (var_ratio,2)

# Store data in a list:
list		= [var_accession, var_ori_reads, var_num_reads, var_dropped,
		var_length, var_GC, var_contigs, var_genus, var_align_genus,
		var_species, var_align_species, var_avg_length, var_ratio, 
		var_coverage, var_date]

print (f"Sequence and assembly analysis:")
print (f"\tAccession:\t\t\t{var_accession}")
print (f"\tOriginal reads:\t\t\t{var_ori_reads}")
print (f"\tProcessed reads:\t\t{var_num_reads}")
print (f"\tDropped reads:\t\t\t{var_dropped} %")
print (f"\tAssembly length:\t\t{var_length} bp")
print (f"\tGC:\t\t\t\t{var_GC} %")
print (f"\tContigs:\t\t\t{var_contigs}")
print (f"\tGenus:\t\t\t\t{var_genus}")
print (f"\tKraken2 genus alignment:\t{var_align_genus} %")
print (f"\tSpecies:\t\t\t{var_species}")
print (f"\tKraken2 species alignment:\t{var_align_species} %")
print (f"\tAvg genome length:\t\t{var_avg_length} bp")
print (f"\tGenome ratio:\t\t\t{var_ratio}")
print (f"\tCoverage:\t\t\t{var_coverage} x")
print (f"\tDate analyzed:\t\t\t{var_date}")

# Print list to file:


# --------------------------------------------------------------------------- #
