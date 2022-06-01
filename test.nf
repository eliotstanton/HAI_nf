#!/usr/bin/env nextflow

// TODO: Add author information

//TODO: Test/fix running only a single pair of FASTQ files
Channel
	.fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
	.set { raw_reads }

//phiX	= 'phiX.fa'
//phiX.view()

python_ch	= Channel.fromPath('data/test.py')

// TODO: Pre-processing step:
// TODO: Ensure that reads are correctly formatted to use HAI-Seq ID nomenclature
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
		ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
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
	tuple name, file("${name}_filtered{_R1,_R2}.fastq") into read_files_spades, read_files_kraken2, read_files_quast

	"""

	ls /data

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

	"""

}

// TODO: Resolve output files not working right
process kraken2 {

	tag "$name"
        publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from read_files_kraken2

	output:
	tuple name, file("${name}_kraken2_report.txt") into files_deliverables

	"""
	echo ${name}
	kraken2 \
		--db /kraken2-db/minikraken2_v1_8GB \
		--threads ${task.cpus} \
		--output ${name}.kraken \
		--use-names \
		--report ${name}_kraken2_report.txt \
		--paired ${name}_filtered_R1.fastq ${name}_filtered_R2.fastq

	"""


}


//TODO: Add producing a trimmed down FASTA file
process spades {

	tag "$name"
        publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(reads) from read_files_spades

	output:
	tuple name, file("${name}.fa") into assembled_genome_mlst, assembled_genome_quast

	"""

	spades.py \
		-1 ${name}_filtered_R1.fastq \
		-2 ${name}_filtered_R2.fastq \
		-o ./output \
		--cov-cutoff 2 \
		--careful \
		--threads ${task.cpus}

	mv ./output/contigs.fasta ${name}.fa

	"""

}

process mlst {

	tag "$name"
        publishDir "${params.outdir}/${name}", mode: 'copy'

	input:
	set val(name), file(assembly) from assembled_genome_mlst

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
	set val(name), file(reads) from read_files_quast
	set val(name2), file(assembly) from assembled_genome_quast

	output:
	tuple name, file("report.txt") into results_quast

	"""
	quast.py \
		--threads ${task.cpus} \
		--output-dir ./quast \
		-1 ${name}_filtered_R1.fastq -2 ${name}_filtered_R2.fastq \
		${name2}.fa \
	"""

}

// TODO: Get and process deliverables
process deliverables {

	input:
//	file test from python_ch

	output:

	"""
	

	"""

}
