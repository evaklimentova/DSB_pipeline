# Pipeline for analysis and visualization of double-strand breaks reparation

This repository contains a pipeline developed for a project dealing with DNA double-strand break reparation.

## Installation guide

1. Clone this repo:

	`git clone https://github.com/evaklimentova/DSB_pipeline.git`

2. Create a conda environment from yml file with all prerequisities:

	`conda env create -n DSB_pipeline -f environment.yml `
	
3. Activate the environment and go to the `source_code` directory:

	`conda activate DSB_pipeline`
	
	`cd source_code`

4. You're prepared to run the pipeline.


## Running the analysis

There are three commands you can run automatically:

- `./generate_genome_index.sh -g <REFERENCE_GENOME_FASTA> -d <DIRECTORY>`

  	This is for generating genome index with STAR, which is then used in the second, mapping, step. As the reference sequence is usually quite short, the script automatically scales down the `genomeSAindexNbases` STAR parameter. 

- `make single_end SAMPLE=<SAMPLE_FASTQ> GENOME_DIR=<GENOME_INDEX_DIR> PREFIX=<PREFIX>`

  	This command will analyse one single end read file. `<SAMPLE_FASTQ>` is the name of your file with reads to analyse, `<GENOME_INDEX_DIR>` is the file, where you have generated the genome indexes from step 1 and `<PREFIX>` prefix name for the output files (default to 'output'). The output files will be stored in the `outputs` directory with the name starting with prefix.

- `make paired_end FORWARD=<FASTQ_F> REVERSE=<FASTQ_R> GENOME_DIR=<GENOME_INDEX_DIR> PREFIX=<PREFIX>`

  	This command is for analysis of two paired end samples. The usage is same as in the analysis of single end sample, but for the input you have to provide 2 FASTQ files - `<FASTQ_F>` for forward reads and `<FASTQ_R>` for reverse reads.
  
 
The `source_code/additional_scripts` directory contains some additional scripts that were later used for an additional analysis.

 

## Examples of running the pipeline

In the `example_inputs` directory, there are prepared some example files for analysis. You can run the following analysis:

- `./generate_genome_index.sh -g ../example_inputs/example_genome-SMG7_amplicon.fas -d ../example_inputs`

  	generating genome index
  
- `make single_end SAMPLE=../example_inputs/example_single_end_reads.fastq GENOME_DIR=../example_inputs PREFIX=example_single_end`

  	analysis of single end data
  
- `make paired_end FORWARD=../example_inputs/example_paired_end_reads-R1.fastq REVERSE=../example_inputs/example_paired_end_reads-R2.fastq GENOME_DIR=../example_inputs PREFIX=example_paired_end`

  	analysis of paired end data


## Outputs

All outputs from all intermediate steps are included, but the most useful and informative files are:

- html report containing summarization of the analysis. It includes coverage plot and viualization of number of insertions and deletions on each position of the reference.

- `<PREFIX>-collapsed.sorted.sam` with aligned reads filtered to contain at least one insertion/deletion and collapsed based on their mapping profile (it's based on the position of deletions and insertions in the reference but it does not care about the read itself; if one read maps to position 1 of the reference and second read maps to position 5 of the reference but both reads have deletion od the same part of the reference, they are collapsed). Each read in the collapsed file is a representative of the whole cluster and it carries a custom `XN` flag showing the original number of reads in this cluster. This file can be viewed in a genome browser such as IGV.






