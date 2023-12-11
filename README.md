# FGF14 analyses

This [Snakemake](https://github.com/snakemake/snakemake) workflow performs the anlysis of long read (Nanopore) repeat sequence data. It calculates statistics and generates plots (custom Python and R scripts) to characterize the nature of the repeat expansion in length and motif composition. Specific alleles (length and motif) can be defined manually for enhanced visualization. 

## Input
The workflow needs fastq.gz files as input. It was tested on fastq files generated from the [FGF14 basecalling](https://github.com/kilpert/FGF14_basecalling) workflow, but this is no requirement.

## Prerequisits
The workflow needs the [Biopyton](https://biopython.org/) library for computing the reverse complement of sequences. Install it beforehand, e.g. in a conda environment, so that it is available when running the workflow.

## Starting the workflow
The workflow can be started by running the `run.sh` on the shell.

## Publication

Mohren et al. (2024). Advancing molecular, phenotypic and mechanistic insights of FGF14 pathogenic expansions (SCA27B). *in prep.*