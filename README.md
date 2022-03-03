# flu-ngs

Influenza virus next generation sequence analysis pipeline.


## Set up

Requirements

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used for generating quality control reports of the raw reads. Developed with 0.11.9. On debian do: `$ apt install fastqc=0.11.9+dfsg-4`.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim adapters off reads. Developed with version 0.39. On debian do `$ apt install trimmomatic=0.39+dfsg-2`.

## Usage

Reads should be placed in a `raw` directory with the following structure. It is fine if the `fastq` files are gzipped (i.e. have a `.gz` suffix).

```
raw/
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

Then run trimming and quality control on these samples:

```bash
snakemake --snakefile workflow/trim-qc.smk -c 8
```

(`-c 8` tells snakemake to use 8 cores, if they're available, so tweak this as you see fit.)

And inspect the HTML output in `results/qc-raw` and `results/qc-trimmed`.