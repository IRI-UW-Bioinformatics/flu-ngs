# flu-ngs

Influenza virus next generation sequence analysis pipeline.


## Requirements

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used for generating quality control reports of the raw reads.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim adapters off reads.
- [IRMA](https://wonder.cdc.gov/amd/flu/irma/) is used to match reads to flu reference sequences.
- [VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) is used to identify effect of nucleotide changes at the protein level.
  - VEP is written in Perl and requires a module called Bundle::DBI. Install it with `perl -MCPAN -e 'install Bundle::DBI'`
- [tabix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/) is required to preprocess files for VEP.
- bgzip & gunzip are used for (de)compression.

Versions are listed in `workflow/envs/*.yaml`.

## Using this repo

Each time you have samples to run, I would suggest cloning this repository:

```bash
git clone git@github.com:IRI-UW-Bioinformatics/flu-ngs.git <name>
```

where `<name>` is the name of the directory that you want, then `cd <name>`.

Reads should be placed in a `raw` directory with the following structure. It is fine if the `fastq` files are gzipped (i.e. have a `.gz` suffix). It is expected that the forward and reverse reads will be in `{sample}_1.fastq` and `{sample}_2.fastq` files. Either rename files accordingly or ask David to program more flexibility. `trimlog.fas` should contain the adapters.

```
raw
├── trimlog.fas
├── YK_2832
│   ├── YK_2832_1.fastq
│   └── YK_2832_2.fastq
├── YK_2833
│   ├── YK_2833_1.fastq
│   └── YK_2833_2.fastq
├── YK_2834
│   ├── YK_2834_1.fastq
│   └── YK_2834_2.fastq
...
```

Specify sample names in a file called `config.json` in the root directory.

```
{
  "samples": [
    "YK_2832",
    "YK_2833",
    "YK_2834"
  ],
  "pair": [
    "combined"
  ]
}
```

`"pair": ["combined"]` tells the pipeline to analyse paired and unpaired reads together.
If you wanted to also run, say, paired reads alone, and unpaired reads alone you would use `"pair": ["combined", "paired", "unpaired"]`.

## Quality control

Run trimming and quality control on these samples:

```bash
snakemake --snakefile workflow/trim-qc.smk --cores all
```

Change the value of `cores` as you like.

Inspect the HTML output in `results/qc-raw` and `results/qc-trimmed`.

## IRMA

Based on the QC results, you could at this point remove samples from `config.json` if you don't want IRMA to waste time analysing them.

Run the IRMA step using:

```bash
snakemake --snakefile workflow/irma.smk --cores all
```

When finished three summary files are generated:

- `results/xlsx/variants-mcc-by-sample-ordered.xlsx`. Each _sample_ has its own sheet. Each sheet contains all flu segments found in that sample.
- `results/xlsx/variants-mcc-by-segment-ordered.xlsx`. Like above, but each _segment_ gets its own sheet.
- `results/xlsx/variants-mcc-flat-ordered.xlsx`. This contains all variants in one flat sheet.

## Sorted BAM files

Most software to look at reads alignments require sorted bam files, and/or bam index files.
I've written a small workflow for generating these for all bam files IRMA generates.
It requires [samtools](http://www.htslib.org/).
`.sorted.bam` and `.sorted.bam.bai` files are saved in the same directory as the original `.bam` files.
Do:

```bash
snakemake --snakefile workflow/sort-bam.smk --cores all
```
