# flu-ngs

Influenza virus next generation sequence analysis pipeline.


## Set up

Requirements

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used for generating quality control reports of the raw reads. On debian do: `$ apt install fastqc=0.11.9+dfsg-4`.

## Usage

Reads should be placed in a `reads` directory with the following structure. It is fine if the `fastq` files are gzipped (i.e. have a `.gz` suffix).

```
reads/
├── trimlog.fas
├── YK_2832
│   ├── YK_2832_1.fastq
│   ├── YK_2832_2.fastq
│   └── YK_2832_final_bwa.bam  # Optional
├── YK_2833
│   ├── YK_2833_1.fastq
│   ├── YK_2833_2.fastq
│   └── YK_2833_final_bwa.bam  # Optional
└── YK_2834
    ├── YK_2834_1.fastq
    ├── YK_2834_2.fastq
    └── YK_2834_final_bwa.bam  # Optional
```

Specify sample names in a file called `config.json`:

```
{
  "samples": [
    "YK_2832",
    "YK_2833",
    "YK_2834"
  ]
}
```

Then run quality control on these samples using:

```bash
snakemake --snakefile workflow/quality-control.smk -c 1
```

And inspect the HTML output in `results/qc`.