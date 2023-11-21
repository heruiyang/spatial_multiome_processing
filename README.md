# spatial_multiome: a Snakemake pipeline to process spatial CUT&Tag and RNA-seq datasets

## Table of Contents

1. [Installation](#installation)
2. [Running the workflow](#running-the-workflow)

## Installation

First, clone the main branch of the repository to your local machine.

```
git clone https://github.com/heruiyang/spatial_multiome_processing/ .
cd spatial_multiome_processing
```

The requirements to run the workflow can be found in `environment.yml`. Create a conda environment with the required dependencies and activate it.

```
conda env create -f environment.yml
conda activate spatial_multiome
```

You're ready to start!

## Setting up a run

Input files and options to the snakemake pipeline are specified in the `config.yaml` config file.

#### Input data

There are two required input `fastq.gz` files and one optional `json` file (work in progress - this is currently required as well, and must be generated from 10X SpaceRanger's manual fiducial alignment). The fastq files are specified in the config file as:

```
sample_fastqs:
  sample_name:
    # Read 1 (16bp Visium barcode + 12bm UMI)
    R1: path/to/read1.fastq.gz
    # Read 2 (gDNA only)
    R2: path/to/read2.fastq.gz
```

The `json` file specifies which spots are part of the tissue and is used to create the filtered feature matrix files in the directory `filtered_peak_bc_matrix`. The pipeline expects an input json file the 10X SpaceRanger format (see [10X Genomics documentation](https://www.10xgenomics.com/support/software/space-ranger/analysis/inputs/image-fiducial-alignment)). This file is specified in the config file as:

```
sample_alignment_jsons:
  # Path to SpaceRanger fiducial alignment json - if you want to do automated alignment (currently not implemented), put 'none' 
  sample_name: path/to/alignment.json
```

#### Reference genome

`spatial_multiome` uses the `bowtie2` aligner and expects references in a compatible format. If you have a `genome.fa` or `genome.fa.gz` file that contains the genome you would like to align to, run the following command to generate a `bowtie2`-compatible reference:

```
bowtie2-build genome.fa {prefix}
```

Replace `prefix` with a suitable name. References are specified in the config file using:

```
# Path and prefix for bowtie2 genome reference
ref: path/to/ref/and/prefix
```

## Running the workflow

Once the `config.yaml` file has been properly configured, run the following command to start the workflow, replacing `{num_cores}` with an integer number of threads to use for execution:

```
snakemake -s workflow.snakefile -j {num_cores}
```

(Optional) before running the full workflow, you can run snakemake with options `-np` to execute a dry run:

```
snakemake -s workflow.snakefile -np
```


