# FGF14 analyses

This [Snakemake](https://github.com/snakemake/snakemake) workflow performs the anlysis of long read (Nanopore) repeat sequence data. It calculates statistics and generates plots (custom Python and R scripts) to characterize the nature of the repeat expansion in length and motif composition. Specific alleles (length and motif) can be defined manually for enhanced visualization. 


## Install from Github

`git clone https://github.com/kilpert/FGF14_analyses.git`


## Setup (Configuration)

The workflow needs **fastq.gz** files as input. It was tested on fastq files generated from the [FGF14 basecalling](https://github.com/kilpert/FGF14_basecalling) workflow, but this is no requirement.

The path to the folder holding the fastq.gz files has to be configured in the `config/config.yaml` file, e.g.:

`fastq_dir: path/to/fastq`


## Starting the workflow
The workflow can be started by running the `run.sh` on the shell 

### or

by executing the Snakemake command directly:

```
snakemake --cores --use-conda --conda-frontend mamba -p --rerun-incomplete
```

All required software packages will be installed via [conda](https://conda.io)/[mamba](https://github.com/mamba-org/mamba) on the first run.


## Demo

A demo dataset (`demo/input/demo.fastq.gz`, 10 reads only) is provided. It can be used for testing the workflow. The run time is about 15s.

Run the workflow on the demo data:

```
snakemake --cores --use-conda --conda-frontend mamba -p --rerun-incomplete --configfile config/demo.config.yaml
```

The demo output will be saved to the `demo/output` folder as specified in the `demo/demo.config.yaml`.

The exact demo output is already provided in the `demo/output.tar.gz` file. Unzip to restore the original output folder from the workflow (e.g. for comparison).


## System requirements
The workflow was tesed on:

 - Ubuntu 22.04.2 LTS
 - conda 23.3.1
 - mamba 1.4.2
 - biopython 1.81
 - Snakemake 7.32.4


## Publication

Mohren et al. (2024). Advancing molecular, phenotypic and mechanistic insights of FGF14 pathogenic expansions (SCA27B). *in prep.*
