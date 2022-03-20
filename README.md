# flu-ngs

Influenza virus next generation sequence analysis pipeline.


## Requirements

### snakemake

The 'pipeline' is really two
[snakemake](https://snakemake.readthedocs.io/en/stable/) workflows:
[`workflow/trim-qc.smk`](workflow/trim-qc.smk), which trims reads and runs
quality control measures, and [`workflow/irma.smk`](workflow/irma.smk) which
runs IRMA and generates summary output.

### Bioinformatics

These workflows call various other bioinformatics programs:

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used
  for generating quality control reports of the raw reads.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim
  adapters off reads.
- [IRMA](https://wonder.cdc.gov/amd/flu/irma/) is used to match reads to flu
  reference sequences.
- [VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) is used to
  identify effect of nucleotide changes at the protein level.
  - VEP is written in Perl and requires a module called Bundle::DBI. Install it
    with `perl -MCPAN -e 'install Bundle::DBI'`
- [tabix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/) is required to
  preprocess files for VEP.
- bgzip & gunzip are used for (de)compression.

Versions are listed in `workflow/envs/*.yaml`.

### Python

There are several python scripts in [`workflow/scripts`](workflow/scripts) which
have a couple of dependencies.

Make a python virtual environment and install them with:

```bash
pip install -r requirements.txt
```

Read more about what python virtual environments are, why they are useful, and
how to set them up
[here](https://realpython.com/python-virtual-environments-a-primer/).

You could use the same virtual environment for each analysis. If you have one
setup, then activate it with:

```bash
source ~/.virtualenvs/flu-ngs-env/bin/activate
```

## Running the pipeline

Each time you have samples to run, I would suggest cloning this repository:

```bash
git clone git@github.com:IRI-UW-Bioinformatics/flu-ngs.git <name>
```

where `<name>` is the name of the directory that you want, then `cd <name>`.

Then, put reads in a directory called `raw` with the following structure:

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

It is fine if the `fastq` files are gzipped (i.e. have a `.gz` suffix).

Forward and reverse reads should be in `{sample}_1.fastq` and `{sample}_2.fastq`
respectively. Either rename files accordingly or ask David to program more
flexibility.

`trimlog.fas` should contain the adapters.

## Quality control

The first step is to trim adapters and run quality control. Specify sample names
in a file called `qc-config.json` in the root directory.

```
{
  "samples": [
    "YK_2832",
    "YK_2833",
    "YK_2834"
  ],
  "pair": [
    "paired",
    "unpaired"
  ]
}
```

`"pair": ["paired", "unpaired"]` tells the pipeline to generate quality control
reports for paired and unpaired reads separately. You could also add
`"combined"` to this array, which would generate one QC report of the paired and
unpaired reads together. This array could also contain only `"combined"`. Do
whatever makes sense for your analysis.

Then, run trimming and quality control on these samples by calling snakemake:

```bash
snakemake --snakefile workflow/trim-qc.smk --cores all
```

Change the value of `cores` as you like.

HTML QC reports are saved in `results/qc-raw` and `results/qc-trimmed`.

## IRMA

A config file `irma-config.json` is required to run IRMA. It could be identical
to `qc-config.json`, but you could remove samples if any didn't pass QC.

If you had set `"pair": ["combined", "paired", "unpaired"]` to look at different
the quality of different types of reads, you may want to now set `"pair":
["combined"]` to just analyse the combined paired and unpaired reads, rather
than conduct three analyses in parallel.

Run the IRMA step using:

```bash
snakemake --snakefile workflow/irma.smk --cores all
```

When finished three summary files are generated:

- `results/xlsx/variants-mcc-by-sample-ordered.xlsx`. Each _sample_ has its own
  sheet. Each sheet contains all flu segments found in that sample.
- `results/xlsx/variants-mcc-by-segment-ordered.xlsx`. Like above, but each
  _segment_ gets its own sheet.
- `results/xlsx/variants-mcc-flat-ordered.xlsx`. This contains all variants in
  one flat sheet.

---

## Bonus: sorted BAM files

Most software to look at reads alignments require sorted bam files, and/or bam
index files. I've written a small workflow for generating these for all bam
files IRMA generates. It requires [samtools](http://www.htslib.org/).
`.sorted.bam` and `.sorted.bam.bai` files are saved in the same directory as the
original `.bam` files. Do:

```bash
snakemake --snakefile workflow/sort-bam.smk --cores all
```
