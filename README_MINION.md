#flu-minion
flu-minion is an Influenza virus next generation sequence analysis pipeline for Oxford MinION data. It is largely similar to flu-ngs except the following changes:

## Requirements

### snakemake 

[`workflow/filter-trim-qc.smk`](workflow/filter-trim-qc.smk) filters and trims reads, 

### Bioinformatics :
- fastqc and trimmomatic steps are not required.
- [chopper] (https://github.com/wdecoster/chopper) is used to filter and trim reads.

### R
There is one R script in [`workflow/scripts`](workflow/scripts) which requires the following packages.
				"data.table", 
                   		"futile.logger",
                   		"ggplot2",
                  		"optparse",
                   		"plyr",
                  		"readr",
                   		"reshape2",
                   		"scales",
                   		"viridis"

## Running the flu-minion workflow

## Combine fastq files
The MinION instrument demultiplexes reads, removes primers and outputs multiple zipped fastq files for each sample. 
The first step is to combine them into a single fastq.gz file. 

The multiple zipped fastq files have to be  put in a directory called `raw` with the following structure:

```
raw
├── barcode05
│   ├── *_0.fastq.gz
│   ├── *_1.fastq.gz
│   ├── *_2.fastq.gz
│   ├── *_3.fastq.gz
│   │...
├── barcode06
│   ├── *_0.fastq.gz
│   ├── *_1.fastq.gz
│   ├── *_2.fastq.gz
│   ├── *_3.fastq.gz
│   │...
├── barcode07
│   ├── *_0.fastq.gz
│   ├── *_1.fastq.gz
│   ├── *_2.fastq.gz
│   ├── *_3.fastq.gz
│   │...
...
```

## filter fastq reads
Then the reads are filter and trim reads based on a quality score and length.
[chopper] (https://github.com/wdecoster/chopper) is used to do this. Currently the minimum quality score is set at 10, minimum read length is 600 and maximum read length is 2500. This could be made flexible for future.

Specify sample names in a file called trim-qc-config.json in the root directory.

```
{
  "samples": [
    "barcode05",
    "barcode06",
    "barcode07"
  ]
}
```
Combine and filter reads in a single step by calling [`workflow/filter-trim-qc.smk`](workflow/filter-trim-qc.smk) 

```bash
snakemake --snakefile workflow/filter-trim-qc.smk --cores all
```

Both the combined file and the filtered file are saved in `trimmed`.

## Minion Quality Control
[Minion Quality Control] (https://github.com/roblanf/minion_qc) generates a range of diagnostic plots and data for quality control of sequencing data [sequencing_summary_*.txt]from Oxford Nanopore's MinION. It is run on [sequencing_summary_*.txt]

call 
```bash
Rscripts workflow/minion-qc.r sequencing_summmary_*.txt
```
