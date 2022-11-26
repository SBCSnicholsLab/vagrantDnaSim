# vagrantDnaSim

Pipeline to simulate vagrant DNA insertion in to nuclear genomes.

## Setup
1. Clone the repository
2. Set up a conda environment: `conda create -n vagrantSim -c bioconda -c conda-forge bwa-mem2 "r-base>=4" r-matrix wgsim samtools freebayes scikit-allel zarr`
3. Run pipeline, for instance: `bash 00pipeline.sh 250M000 250000000 0.00004`. The parameters are: [outname] [nuclear genome size] [insertion rate per year]

## What it does

1. In R
  * Generates an empty nuclear genome of defined size (sparse array to save memory).
  * Genetates a random extranuclear sequence of 16000nt length
  * Cylces through 15,000,000 years, checking each year if a substitution happens in the extranuclear DNA and is an insertion happens into the nuclear genome. 
  * Genomes are written out. The nuclear one is made up randomly while on the fly, filling in vagrant inserts that were recorded during the 15M years.
2. Sequencing reads are simulated using [wgsim](https://github.com/lh3/wgsim)
