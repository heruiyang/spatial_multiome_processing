import yaml
import os
import pandas as pd
import re

####	Set up config

configfile: "config.yaml"

input_samples=config["sample_fastqs"]


####	End of config


### Process

rule all:
	input:
		expand("results/{sample}/fastq_process/{sample}_processed.barcoded.fastq.gz", sample=input_samples.keys()),
		expand("results/{sample}/{sample}_frags.sorted.bed.gz", sample=input_samples.keys()),
		expand("results/{sample}/{sample}_frags.sorted.bed.gz.tbi", sample=input_samples.keys()),
		expand("results/{sample}/align/{sample}_bc.sorted.dedup.bam", sample=input_samples.keys()),
		expand("results/{sample}/align/{sample}_bc.sorted.dedup.bam.bai", sample=input_samples.keys()),
		expand("results/{sample}/macs2_callpeak/{sample}_peaks.narrowPeak", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}.sorted.dedup.bam", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_edit_distance.tsv", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_per_umi_per_position.tsv", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_per_umi.tsv", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/fastqc/{sample}_r1_trimmed_fastqc.html", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/fastqc/{sample}_r2_trimmed_fastqc.html", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/cutadapt/{sample}_r1_discarded.fastq.gz", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/cutadapt/{sample}_r2_discarded.fastq.gz", sample=input_samples.keys())


'''rule process_bc_samples:
	input:
		read1_fastq = expand(input_r1)
	output:
		umi_fastq="results/{sample}_umis.fastq",
		bc_fastq="results/{sample}_barcodes.fastq"
	message:
		"Splitting read 1 fastqs"
	shell:
		"zcat {input.read1_fastq} | python scripts/split_fastq.py - {output.umi_fastq} {output.bc_fastq}"
'''

rule fastq_trim:
	input:
		read1 = lambda wildcards: config['sample_fastqs'][wildcards.sample]['R1'],
		read2 = lambda wildcards: config['sample_fastqs'][wildcards.sample]['R2']
	output:
		read1_trimmed = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2_trimmed = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz",
		read1_discarded = "results/{sample}/fastq_process/cutadapt/{sample}_r1_discarded.fastq.gz",
		read2_discarded = "results/{sample}/fastq_process/cutadapt/{sample}_r2_discarded.fastq.gz",
		discarded_read_names = "results/{sample}/fastq_process/cutadapt/discarded_reads.lst"
	message:
		"Trimming R1 and R2"
	log:
		"logs/{sample}/fastq_trim.log"
	threads: 8
	run:
		shell("cutadapt --minimum-length 28 --cores={threads} -a 'AAAAAAX;o=6' -a 'TTTTTTX;o=6' -a 'CCCCCCX;o=6' -a 'GGGGGGX;o=6' -A 'AAAAAAX;o=6' -A 'TTTTTTX;o=6' -A 'CCCCCCX;o=6' -A 'GGGGGGX;o=6' --times 15 -o {output.read1_trimmed} -p {output.read2_trimmed} {input.read1} {input.read2} &> {log}")
		shell("comm -23 <( zcat {input.read1} | paste - - - - | grep -Po '^@\K\S+' | sort ) <( zcat {output.read1_trimmed} | paste - - - - | grep -Po '^@\K\S+' | sort ) > {output.discarded_read_names}")
		shell("seqtk subseq {input.read1} {output.discarded_read_names} | gzip > {output.read1_discarded}")
		shell("seqtk subseq {input.read2} {output.discarded_read_names} | gzip > {output.read2_discarded}")

rule fastqc:
	input:
		read1_trimmed = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2_trimmed = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz"
	output:
		read1_report = "results/{sample}/fastq_process/fastqc/{sample}_r1_trimmed_fastqc.html",
		read2_report = "results/{sample}/fastq_process/fastqc/{sample}_r2_trimmed_fastqc.html"
	message:
		"Generating FastQC read quality report"
	log:
		"logs/{sample}/fastqc.log"
	run:
		shell("fastqc -o results/{wildcards.sample}/fastq_process/fastqc -f fastq {input.read1_trimmed} {input.read2_trimmed} &> {log}")

rule process_umi_bc:
	input:
		read1 = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2 = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz"
	output:
		umitools_fastq=temporary("results/{sample}/fastq_process/{sample}_processed.fastq.gz"),
		fastq="results/{sample}/fastq_process/{sample}_processed.barcoded.fastq.gz"
	message:
		"Filtering barcodes and performing UMI deduplication"
	params:
		whitelist=config["whitelist"]
	log: 
		umitools_log="logs/{sample}/process_umi_bc_umitools.log",
		sinto_log="logs/{sample}/process_umi_bc_sinto.log"
	run:
		shell("zcat {input.read1} | umi_tools extract --extract-method=string --bc-pattern=XXXXXXXXXXXXXXXXNNNNNNNNNNNN --read2-in={input.read2} --read2-out={output.umitools_fastq} --log {log.umitools_log} 1>/dev/null")
		shell("sinto barcode --barcode_fastq {input.read1} --read1 {output.umitools_fastq} -b 16 --whitelist {params.whitelist} &>{log.sinto_log}")


rule align:
	input:
		fastq = "results/{sample}/fastq_process/{sample}_processed.barcoded.fastq.gz"
	output:
		bam = "results/{sample}/align/{sample}.sorted.bam",
		index = "results/{sample}/align/{sample}.sorted.bam.bai"
	message:
		"Aligning to specified reference genome with bowtie2"
	params:
		ref=config["ref"]
	log: "logs/{sample}/align.log"
	threads: 8
	run:
		shell("bowtie2 --very-sensitive -q --phred33 --end-to-end -t -x {params.ref} -U {input.fastq} -p {threads} | samtools view -bS - | samtools sort - -o {output.bam} -@ {threads} &>{log}")
		shell("samtools index {output.bam} -o {output.index} &>{log}")


rule dedup:
	input:
		bam = "results/{sample}/align/{sample}.sorted.bam"
	output:
#		statsfile1 = "results/{sample}/{sample}_edit_distance.tsv",
#		statsfile2 = "results/{sample}/{sample}_per_umi_per_position.tsv",
#		statsfile3 = "results/{sample}/{sample}_per_umi.tsv",
		dedup_bam = temporary("results/{sample}/align/{sample}.sorted.dedup.bam"),
		bam = "results/{sample}/align/{sample}_bc.sorted.dedup.bam",
		index = "results/{sample}/align/{sample}_bc.sorted.dedup.bam.bai"
	message:
		"Performing UMI deduplication with umi_tools and bam file barcode assignment"
	log: 
		umitools_log = "logs/{sample}/dedup_umitools.log",
		sinto_log = "logs/{sample}/dedup_sinto.log"
	params:
		use_umi = config['umi_dedup']
	run:
		if params.use_umi:
			shell("umi_tools dedup --stdin={input.bam} --log={log.umitools_log} --per-cell --output-stats=results/{wildcards.sample}/{wildcards.sample} > {output.dedup_bam}")
		else:
			shell("umi_tools dedup --stdin={input.bam} --log={log.umitools_log} --per-cell --ignore-umi > {output.dedup_bam}")
		shell("samtools view -h {output.dedup_bam} | sinto nametotag -b - | samtools view -b - > {output.bam} 2>{log.sinto_log}")
		shell("samtools index {output.bam} -o {output.index}")


"""
Note: our data is single end so we fake a fragments file by just converting our bam file to a correctly formatted tsv
according to https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments

We should also consider using e.g. MACS2 estimated fragment length to generate a more correct fragments file
"""
rule build_fragments_file:
	input:
		bam = "results/{sample}/align/{sample}_bc.sorted.dedup.bam"
	output:
		frag_file_unsorted = temporary("results/{sample}/{sample}_frags.bed"),
		frag_file_sorted = temporary("results/{sample}/{sample}_frags.sorted.bed"),
		frag_file_sorted_gzip = "results/{sample}/{sample}_frags.sorted.bed.gz",
		frag_file_sorted_index = "results/{sample}/{sample}_frags.sorted.bed.gz.tbi"
	log: "logs/{sample}/build_fragments_file.log"
	message: "Creating fragments file"
	run:
		shell("python ./scripts/create_fragments.py {input.bam} {output.frag_file_unsorted} 2>{log}")
		shell("sort -k1,1 -k2,2n {output.frag_file_unsorted} 1>{output.frag_file_sorted} 2>{log}")
		shell("bgzip {output.frag_file_sorted} -c 1>{output.frag_file_sorted_gzip} 2>{log}")
		shell("tabix -p bed {output.frag_file_sorted_gzip} 2>{log}")


rule call_peaks:
	input:
		bam = "results/{sample}/align/{sample}_bc.sorted.dedup.bam"
	output:
		peaks = "results/{sample}/macs2_callpeak/{sample}_peaks.narrowPeak"
	message: "Calling peaks with MACS2"
	log: "logs/{sample}/call_peaks.log"
	params:
		outdir = "results/{sample}/macs2_callpeak",
		fdr = config['macs2_fdr'],
		genome_size = config['macs2_genomesize'],
		name = '{sample}'
	shell:
		"""	
		macs2 callpeak --SPMR -B -q {params.fdr} --keep-dup all \
           -g {params.genome_size} -f BAM \
           -t {input.bam} --outdir {params.outdir} -n {params.name} > {log} 2>&1
        """


rule create_feature_mtx:
	input:
		frag_file = "results/{sample}/{sample}_frags.sorted.bed.gz",
		peaks = "results/{sample}/macs2_callpeak/{sample}_peaks.narrowPeak"
	output:
		"results/{sample}/raw_peak_bc_matrix/barcodes.tsv",
		"results/{sample}/raw_peak_bc_matrix/peaks.bed",
		"results/{sample}/raw_peak_bc_matrix/matrix.mtx"
	message: "Creating feature by spot matrix"
	log: "logs/{sample}/create_feature_mtx.log"
	params:
		output_dir = "results/{sample}/raw_peak_bc_matrix"
	shell:
		"Rscript ./scripts/create_feature_mtx.R {input.frag_file} {input.peaks} {params.output_dir} &> {log}"


rule filter_feature_mtx:
	input:
		barcodes = "results/{sample}/raw_peak_bc_matrix/barcodes.tsv",
		peaks = "results/{sample}/raw_peak_bc_matrix/peaks.bed",
		matrix = "results/{sample}/raw_peak_bc_matrix/matrix.mtx",
	output:
		barcodes = "results/{sample}/filtered_peak_bc_matrix/barcodes.tsv",
		peaks = "results/{sample}/filtered_peak_bc_matrix/peaks.bed",
		matrix = "results/{sample}/filtered_peak_bc_matrix/matrix.mtx"
	message: "Filtering spots by aligned tissue image"
	log: "logs/{sample}/filter_feature_mtx.log"
	params:
		json = lambda wildcards: config["sample_alignment_jsons"][wildcards.sample],
		output_dir="results/{sample}/filtered_peak_bc_matrix",
		spot_coords=config["spot_coords"]
	run:
		if (params.json != 'none'):
			shell('Rscript ./scripts/filter_feature_mtx.R {params.json} {input.barcodes} {input.peaks} {input.matrix} {params.spot_coords} {params.output_dir} &> {log}')
		else:
			shell('python -u ./scripts/align_image.py')
			shell('Rscript ./scripts/filter_feature_mtx.R results/{sample}/fiducial.json {input.barcodes} {input.peaks} {input.matrix} {params.spot_coords} {params.output_dir} &> {log}')


'''
rule process_bam:
	input:
		bam = "results/{sample}/{sample}_dedup.bam"
	output:
		bam = "results/{sample}/{sample}_dedup_processed.bam"
		index = "results/{sample}/{sample}_dedup_processed.bam.bai"
	message:
		"Processing deduplicated BAM file"
	log: "logs/{sample}/{sample}_process_bam.log"
	shell:
		"samtools index {output.bam} -o {output.index}"
'''


