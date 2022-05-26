#!/usr/bin/env nextflow

// println "Hello, World!"

string_hello	= "Hello, World!"

map_hello	= ["string_hello":15, 1:"Hello", "Hello":100, "World":12, "!":"World"]

var_random	= Math.random()

//println "$string_hello\n"

//println map_hello.Hello

//var_random *= 100

//println var_random

Channel
	.fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
	.into { raw_reads }

phiX	= Channel.fromPath( '/panfs/roc/groups/7/mdh/estanton/HAI_nf/phiX.fa' )
//phiX.view()

// TODO: Pre-processing step:
process preProcess {

	input:
	set val(name), file(reads) from raw_reads

	output:
	tuple name, file(outfiles) into read_files_trimming

	script:
	outfiles	= reads
//	println reads

	"""
	"""

}

// TODO: Job trimming and filtering steps:
process trimming {

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
		ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
		LEADING:20 TRAILING:20 MINLEN:50 AVGQUAL:30

	"""

}

process filtering {

	input:
	set val(name), file(reads) from read_files_filtering
	path '/panfs/roc/groups/7/mdh/estanton/HAI_nf/phiX.fa'

	output:
//	tuple name, file("${name}_filtered{_1,_2}.fastq.gz") into read_files_spades
//	file("${name}_filtered{_1,_2}.fastq.gz") into read_files_kraken2

	"""

	bbmap.sh \
		threads=${task.cpus} \
		ref='/panfs/roc/groups/7/mdh/estanton/HAI_nf/phiX.fa' nodisk \
		in=${name}_trimmed_R1.fastq \
		in2=${name}_trimmed_R2.fastq \
		outu=${name}_filtered.fastq \
		outm=${name}_phiX_R1.fastq

        reformat.sh \
                overwrite=true \
                in=${name}_filtered.fastq \
                out1=${name}_filtered_R1.fastq \
                out2=${name}_filtered_R2.fastq
	"""

}

// Run Kraken2 and Spades concurrently

// Kraken2
//process kraken2 {


//}

// Spades assembly and  MLST
//process spades {


//}

// TODO: Get deliverables
//process deliverables {


//}
