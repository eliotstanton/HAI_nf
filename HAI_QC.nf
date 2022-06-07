#!/usr/bin/env nextflow

// ----------------------------------------------------------------------------

// File:		HAI_QC2.nf
// Date created:	07 June 2022
// Date last modified:	07 June 2022
// Author:		Eliot Stanton (eliot.stanton@state.mn.us)
// Description:		Perform QC analysis, genome assembly, and species ID
//			on WGS data for upload to NCBI.

// ----------------------------------------------------------------------------

// Channels:

Channel
	.fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fa}.gz", size: 2)
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
	.set { raw_reads }

process preProcess {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from raw_reads

	output:
	tuple name, file(outfiles) into read_files_trimming
	tuple name, file(reads) into read_files_deliverables
	
	script:
	outfiles	= reads

	"""

	md5sum ${reads[0]} > ${reads[0]}.md5
	md5sum ${reads[1]} > ${reads[1]}.md5

	"""

}

process trimming {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from read_files_trimming

	output:
	tuple name, file("${name}_trimmed{_R1,_R2}.fastq") into read_files_filtering

	"""

	trimmomatic PE \
		-threads ${task.cpus} \
		-phred33 \
		${reads[0]} ${reads[1]} \
		${name}_trimmed_R1.fastq ${name}_trimmed_unpaired_R1.fastq \
		${name}_trimmed_R2.fastq ${name}_trimmed_unpaired_R2.fastq \
		ILLUMINACLIP:${params.clip} \
		LEADING:${params.leading} \
		TRAILING:${params.trailing} \
		MINLEN:${params.minlen} \
		AVGQUAL:${params.avgqual}

	"""

}

process filtering {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from read_files_filtering

	output:
	tuple name, file("${name}_filtered{_R1,_R2}.fastq.gz") into filtered_spades, 
	filtered_kraken2, filtered_quast, filtered_deliverables

	
	"""

	bbmap.sh \
		threads=${task.cpus} \
		ref='/bbmap/resources/phix174_ill.ref.fa.gz' nodisk \
		in=${name}_trimmed_R1.fastq \
		in2=${name}_trimmed_R2.fastq \
		outu=${name}_filtered.fastq \
		outm=${name}_phiX_R1.fastq

	reformat.sh \
		overwrite=true \
		in=${name}_filtered.fastq \
		out1=${name}_filtered_R1.fastq \
		out2=${name}_filtered_R2.fastq

	gzip ${name}_filtered_R1.fastq
	gzip ${name}_filtered_R2.fastq

	"""
}

process kraken2 {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from filtered_kraken2

	output:
	tuple name, file("${name}_kraken2_report.txt") into kraken2_results

	"""
        
	kraken2 \
		--db /kraken2-db/minikraken2_v1_8GB \
		--threads ${task.cpus} \
		--output ${name}.kraken \
		--use-names \
		--report ${name}_kraken2_report.txt \
		--paired ${name}_filtered_R1.fastq.gz ${name}_filtered_R2.fastq.gz

	"""


}

process spades {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from filtered_spades

	output:
	tuple name, file("${name}.fa") into assembly_mlst, assembly_quast, assembly_contigs

        """

	spades.py \
		-1 ${name}_filtered_R1.fastq.gz \
		-2 ${name}_filtered_R2.fastq.gz \
		-o ./output \
		--cov-cutoff 2 \
		--careful \
		--threads ${task.cpus}

	mv ./output/contigs.fasta ${name}.fa

"""

}

process contigs {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from assembly_contigs

	output:
	tuple name, file("${name}_trimmed.fa")

//	shell:

	"""

	#!/usr/bin/env perl
	use strict;
	use warnings;

	my \$var_minlen = 500;
	my \$file_in    = "${name}.fa";
	my \$file_out   = "${name}_trimmed.fa";

	# Open \$file_in:
	open ( my \$file_read, '<', \$file_in ) or die "Unable to open \$file_in!";

	# Open \$file_out:
	open ( my \$file_write, '>', \$file_out ) or die "Unable to open \$file_out!";

	while ( <\$file_read> ){

		if ( \$_ =~ /^>/ ) {

			my @array_temp  = split ("_", \$_);

			my \$var_length = \$array_temp[3];

			if ( \$var_length < \$var_minlen ) { last }

			else { print \$file_write \$_; }

		}

		else { print \$file_write \$_; }

	}

	# Close \$file_in and \$file_out:
	close ( \$file_read ) && close ( \$file_write );

	"""

}

process mlst {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(assembly) from assembly_mlst

	output:
	tuple name, file("${name}_mlst.txt") into results_mlst

	"""

	mlst \
		--threads ${task.cpus} \
		${name}.fa > ${name}_mlst.txt

	"""

}

process quast {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from filtered_quast
	set val(name2), file(assembly) from assembly_quast

	output:
	tuple name, file("./quast/report.txt") into quast_results

	"""

	quast.py \
		--threads ${task.cpus} \
		--output-dir ./quast \
		-1 ${name}_filtered_R1.fastq.gz -2 ${name}_filtered_R2.fastq.gz \
		${name2}.fa \

	"""

}

process deliverables {

	tag "$name"
	publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from read_files_deliverables
	set val(name2), file(results) from kraken2_results
	set val(name3), file(assembly) from quast_results
	set val(name4), file(reads) from filtered_deliverables

	output:
	tuple name, file("${name}_stats.txt") into deliverables_results

	"""

	#!/usr/bin/env python3

	import sys
	import subprocess
	import pandas as pd
	import os
	import re

	# Variables imported by this script:
	var_accession	= "${name}"

	# Filenames used by this script:
	file_kraken2	= "${name2}_kraken2_report.txt"
	file_FASTQ1	= "${name4}_filtered_R1.fastq"
	file_FASTQ2	= "${name}_R1.fastq.gz"
	file_quast	= "report.txt"
	file_out	= "${name}_stats.txt"

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

	# Determine the number of original paired reads:
	var_ori_reads	= "echo \$(zcat {} | wc -l) /4 | bc".format(file_FASTQ2)
	var_ori_reads	= subprocess.getoutput(var_ori_reads)

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
	data_list_df	= pd.DataFrame(data_list, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
	data_genus_df	= pd.DataFrame(data_genus, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
	data_species_df	= pd.DataFrame(data_species, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])

	# Remove left leading spaces from the Name column
	data_list_df['Name']    = data_list_df['Name'].str.lstrip()
	data_genus_df['Name']   = data_genus_df['Name'].str.lstrip()
	data_species_df['Name'] = data_species_df['Name'].str.lstrip()

	# Sort dataframes by percentages:
	data_list_df	= data_list_df.sort_values(by=['Percentage'], ascending=False)
	data_genus_df	= data_genus_df.sort_values(by=['Percentage'], ascending=False)
	data_species_df	= data_species_df.sort_values(by=['Percentage'], ascending=False)

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
	var_genus	= genus_df['Name'].iloc[0]
	var_species	= species_df['Name'].iloc[0]

	list_perc		= list_df['Percentage'].iloc[0]
	var_align_genus		= genus_df['Percentage'].iloc[0]
	var_align_species	= species_df['Percentage'].iloc[0]
	var_species		= var_species.split(' ')
	var_species		= var_species[-1]

	# Scrape average species genome length from NCBI:
	var_html	= var_html + var_genus + "+" + var_species
	file_html	= var_accession + ".html"
	wget		= "wget --output-document {} -q {}".format(file_html,var_html)
	os.system(wget)

	with open(file_html) as file_read:
		for line in file_read:
			search = re.search("median total length",line)
			if search:
				sline = re.split('\\>|\\<',line)
				for i in sline:
					search = re.search("median total length",i)
					if search:
						j = re.split(': ', i)
						var_avg_length = j[-1]
						var_avg_length = float(var_avg_length)
						var_avg_length *= 1000000
						var_avg_length = int(var_avg_length)

	# Iterate through Quast results:
	with open(file_quast) as file_read:
		for line in file_read:
			# Total assembly length:
			search_total_length = re.search("^Total length \\(>= 0 bp\\)",line)
			if search_total_length:
				var_length	= re.split(' ', line)
				var_length	= var_length[10]
		
			# Number of contigs:
			search_contigs = re.search("# contigs  ", line)
			if search_contigs:
				var_contigs	= re.split (' ', line)
				var_contigs	= var_contigs[20]

			# GC content:
			search_GC = re.search("GC ", line)
			if search_GC:
				var_GC	= re.split (' ', line)
				var_GC	= var_GC[23]

			# Post-processed read number:
			search_left = re.search("# left", line)
			if search_left:
				var_num_reads	= re.split (' ', line)
				var_num_reads	= var_num_reads[23]

			# Read depth:
			search_coverage = re.search("Avg. coverage depth", line)
			if search_coverage:
				var_coverage	= re.split(' ', line)
				var_coverage	= var_coverage[11]

	# Calculate dropped read percentage and round to two decimals:
	var_num_reads   = int(var_num_reads)
	var_ori_reads   = int(var_ori_reads)
	var_dropped     = 1-(var_num_reads/var_ori_reads)
	var_dropped     = var_dropped * 100
	var_dropped     = round(var_dropped, 2)

	# Determine genome ratio and round to two decimals:
	var_length      = int(var_length)
	var_avg_length  = int(var_avg_length)
	if ( var_avg_length >0 ):
        	var_ratio       = var_length/var_avg_length
	var_ratio       = round (var_ratio,2)

	list0	= ['Accession', 'Original reads', 'Processed reads', 'Dropped (%)',
		'Assembled genome length', 'GC (%)', 'Contigs', 'Genus',
		'Kraken2 genus alignment', 'Species', 'Kraken2 species alignment',
		'Avg genome length', 'Genome ratio', 'Coverage depth', 'Date analysed']

	# Store data in a list:
	list1	= [var_accession, var_ori_reads, var_num_reads, var_dropped,
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

	str0    = ""
	str1    = ""

	for element in list0:
	        element = str(element)
	        str0 = str0 + "\t" + element

	for element in list1:
	        element = str(element)
	        str1 = str1 + "\t" + element

	with open(file_out, 'w') as f:
	        f.write(str0)
	        f.write(str1)

	"""

}
