# flu-ngs

Influenza virus next generation sequence analysis pipeline.


## Requirements

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used for generating quality control reports of the raw reads.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim adapters off reads.
- [IRMA](https://wonder.cdc.gov/amd/flu/irma/) is used to match reads to flu reference sequences.
- [VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) is used to identify effect of nucleotide changes at the protein level.
- [tabix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/) is required to preprocess files for VEP.
- bgzip & gunzip are used for (de)compression.

Versions are listed in `workflow/envs/*.yaml`.

## Quality control

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

When finished, a summary of the variants found is saved in an excel sheet in `results/variants/{sample}_{pair}/{sample}_{pair}.xlsx"` where `{pair}` is either `combined` or `paired`.

By default all steps are run for paired and combined (paired and unpaired) reads. This might approximately double the runtime compared to only running, say, combined. This would be fairly trivial to alter.
