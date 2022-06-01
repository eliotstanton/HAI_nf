#!/usr/bin/python3

import sys
import pandas as pd

# --------------------------------------------------------------------------- #

# File name:		test2.py
# Date created:		01 June 2022
# Date last modified:	01 June 2022
# Author:		Eliot Stanton (eliot.stanton@state.mn.us)
# Description:		This script focuses solely on algorithm for pulling
#			species data from Kraken2 output.

# --------------------------------------------------------------------------- #

# Variables imported by this script:
# TODO: Make input of command-line variables fail gracefully
var_accession   = sys.argv[1]
dir_in          = sys.argv[2]
dir_out         = sys.argv[3]

print(f"{var_accession}\n{dir_in}\n{dir_out}\n")

# Define help string:
var_help	= "test2.py [ACCESSION] [DIR_IN] [DIR_OUT]\n"

# Handle missing variables from the command line:
#if ( ):
#	print(var_help)

# Filepaths used by this script:
file_kraken2	= "{}/{}/{}_kraken2_report.txt".format(dir_out, var_accession, var_accession)

# Data structures used by this script:

# --------------------------------------------------------------------------- #

print(f"{var_accession}\n{dir_in}\n{dir_out}\n\n{file_kraken2}")

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
#			print (split_line)

		# Get species results (denoted by 'S') and append to data
		if split_line[3] == 'S':
			data_list.append(split_line)
			data_species.append(split_line)
#			print (split_line)

		# Get species results (denoted by 'G') and append to data
		if split_line[3] == 'G':
			data_list.append(split_line)
			data_genus.append(split_line)
#			print (split_line)


# Convert data lists to data frames: 
data_list_df	= pd.DataFrame(data_list, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
data_genus_df	= pd.DataFrame(data_genus, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
data_species_df	= pd.DataFrame(data_species, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])

# Remove left leading spaces from the Name column
data_list_df['Name']	= data_list_df['Name'].str.lstrip()
data_genus_df['Name']	= data_genus_df['Name'].str.lstrip()
data_species_df['Name']	= data_species_df['Name'].str.lstrip()

# Sort dataframes by percentages:
data_list_df	= data_list_df.sort_values(by=['Percentage'], ascending=False)
data_genus_df	= data_genus_df.sort_values(by=['Percentage'], ascending=False)
data_species_df = data_species_df.sort_values(by=['Percentage'], ascending=False)


# Confine to Name and Percentage columns:
list_df		= data_list_df[['Name', 'Percentage']]
genus_df	= data_genus_df[['Name', 'Percentage']]
species_df	= data_species_df[['Name', 'Percentage']]

# Get first match from each data frame:
list_df		= list_df.head(1) 
genus_df	= genus_df.head(1)
species_df	= species_df.head(1)

# Pull values:
list_name	= list_df['Name'].iloc[0]
genus_name	= genus_df['Name'].iloc[0]
species_name	= species_df['Name'].iloc[0]

list_perc       = list_df['Percentage'].iloc[0]
genus_perc      = genus_df['Percentage'].iloc[0]
species_perc    = species_df['Percentage'].iloc[0]

#genus_df	= data_genus_df['Name']
#genus_df	= genus_df.head(1)



print(f"{list_name}\n{genus_name}\n{species_name}\n")

print(f"{list_perc}\n{genus_perc}\n{species_perc}\n")

species_name	= species_name.split(' ')

species_name	= species_name[-1]

print (genus_name, genus_perc, species_name, species_perc)
