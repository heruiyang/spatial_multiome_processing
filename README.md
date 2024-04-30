# spatial_multiome: a Snakemake pipeline to process spatial CUT&Tag and RNA-seq datasets

## Table of Contents

1. [Installation](#installation)
2. [Setting up a run](#setting-up-a-run)
3. [Configuring workflow options](#configuring-workflow-options)
4. [Running the workflow](#running-the-workflow)

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

#### `sample_fastqs` - Input data

```
sample_fastqs:
  sample_name:
    # Read 1 (16bp Visium barcode + 12bm UMI)
    R1: path/to/read1.fastq.gz
    # Read 2 (gDNA only)
    R2: path/to/read2.fastq.gz
```

#### `sample_alignment_jsons` - Tissue image alignment

```
sample_alignment_jsons:
  # Path to SpaceRanger fiducial alignment json - if you want to do automated alignment (currently not implemented), put 'none' 
  sample_name: path/to/alignment.json
```

The `json` file specifies which spots are part of the tissue and is used to create the filtered feature matrix files in the directory `filtered_peak_bc_matrix`. The pipeline expects an input json file the 10X SpaceRanger format (see [10X Genomics documentation](https://www.10xgenomics.com/support/software/space-ranger/analysis/inputs/image-fiducial-alignment)). 

#### `ref` - Reference genome

`spatial_multiome` uses the `bowtie2` aligner and expects references in a compatible format. If you have a `genome.fa` or `genome.fa.gz` file that contains the genome you would like to align to, run the following command to generate a `bowtie2`-compatible reference:

```
bowtie2-build genome.fa {prefix}
```

Replace `prefix` with a suitable name. References are specified in the config file using:

```
# Path and prefix for bowtie2 genome reference
ref: path/to/ref/and/prefix
```

## Configuring workflow options

The workflow has a number of options that can be adjusted. These options are specified in the `config.yaml` configuration file, where each option is specified as a `param:value` pair. Options are explained below:

#### `whitelist`

```
# Path for text file containing barcode whitelist
whitelist: refs/visium-v1.txt
```

The whitelist text file contains all valid spatial barcodes, with one barcode on each line. For the Visium v1 protocol, this file is provided in this repository.

#### `spot_coords`

```
# Path for text file containing spatial row/column indices of valid barcodes
spot_coords: refs/visium-v1_coordinates.txt
```

The spot coordinates file is a tab-delimited file containing:

Column 1: all valid spatial barcodes (these must match the barcodes specified in `whitelist` above) \
Column 2: 1-based column index of the spatial barcode  \
Column 3: 1-based row index of the spatial barcode 

#### `umi-dedup`

```
# Whether to use UMIs for deduplication
umi_dedup: False
```

Boolean value that specifies whether deduplication should be performed by UMI or by genomic alignment position.

### Peak calling options

#### `macs2_fdr`

```
### Options for MACS2 peak calling
# FDR threshold for peak calling
macs2_fdr: 0.01
```

False discovery rate for MACS2 peak calling. Reasonable options are `0.01`, `0.05`, etc.

#### `macs2_genomesize`

```
# Effective genome size
# Defaults for species:
#   hs: 2.7e9
#   mm: 1.87e9
#   ce: 9e7
#   dm: 1.2e8
macs2_genomesize: 1.87e9
```

Effective genome size for MACS2 peak calling. Default sizes for common species are listed above.

## Running the workflow

Once the `config.yaml` file has been properly configured, run the following command to start the workflow, replacing `{num_cores}` with an integer number of threads to use for execution:

```
snakemake -s workflow.snakefile -j {num_cores}
```

Alternatively, if you are running this in a SLURM-managed HPC system, you can use the `run_workflow_sbatch.sh` script to run the pipeline in the background using sbatch. In that case, run:

```
sbatch run_workflow_sbatch.sh
```

This should work assuming you have set up the environment correctly; if not, check the sbatch parameters at the top of the script.

(Optional) before running the full workflow, you can run snakemake with options `-np` to execute a dry run:

```
snakemake -s workflow.snakefile -np
```


