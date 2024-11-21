# flu-ngs

Influenza virus next generation sequence workflows.

- Uses [IRMA](https://wonder.cdc.gov/amd/flu/irma/) to iteratively refine
  reference sequences.
- Analyses variants for standard influenza A and B segments, and influenza A
  splice variants: PA-X, PB1-F2, M1, M2.
- Reports coding effects of nucleotide changes, taking account of multiple SNPs
  in a codon, and their phase.
- Support for MinION and MiSeq data.
- Quasispecies and mixed sample sequence reconstruction.

## Requirements

### snakemake

The 'pipeline' is really two
[snakemake](https://snakemake.readthedocs.io/en/stable/) workflows:
[`preprocess.smk`](workflow/preprocess.smk), which trims reads and runs
quality control measures, and [`irma.smk`](workflow/irma.smk) which
runs IRMA and generates summary output.

These workflows call various other programs:

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used
  for generating quality control reports of the raw reads.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim
  adapters off reads.
- [IRMA](https://wonder.cdc.gov/amd/flu/irma/) is used to match reads to flu
  reference sequences.
- [VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) is used to
  identify effect of nucleotide changes at the protein level.
  - VEP is written in Perl and requires a module called `Bundle::DBI`. Install it
    with `perl -MCPAN -e 'install Bundle::DBI'`
- [tabix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/) is required to
  preprocess files for VEP. Install via `perl -MCPAN -e 'install Bio::DB::HTS::Tabix'`.
- bgzip & gunzip are used for (de)compression.
- [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml) is used to write
  transcripts from fasta and GFF files.
- [transeq](https://www.ebi.ac.uk/Tools/emboss/) is used to translate
  transcripts.
- [chopper](https://github.com/wdecoster/chopper) is used for filtering and trimming MinION data.
- [MinIONQC](https://github.com/roblanf/minion_qc) is used for generating quality control reports of
  MinION data.
- [TenSQR](https://github.com/davipatti/TenSQR) and [aBayesQR](https://github.com/SoYeonA/aBayesQR/)
  are used for quasispecies / mixed sample sequence inference.
- [minimap2](https://github.com/lh3/minimap2) is used to map reads to IRMA reference sequences for
  quasispecies / mixed sample reconstruction.

Versions are listed in `workflow/envs/*.yaml`.

### Python

There are several python scripts in [`workflow/scripts`](workflow/scripts) which
have a couple of dependencies.

Make a python virtual environment and install the requirements with:

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

## Running the workflow

Each time you have samples to run, I would suggest cloning this repository:

```bash
git clone git@github.com:IRI-UW-Bioinformatics/flu-ngs.git <name>
```

where `<name>` is the name of the directory that you want, then `cd <name>`.

## Data

The next step is to put the read files in a structure expected by the workflow.

### MiSeq

Put reads in a directory called `raw` with the following structure:

```
raw/
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
respectively.

`trimlog.fas` should contain the adapters.

### MinION

Put reads in a directory called `raw`. They must be gzipped (end with `fastq.gz`).

```
raw/
├── barcode05
│   ├── FAW31148_pass_barcode05_485f6488_e81f1340_820.fastq.gz
│   ├── FAW31148_pass_barcode05_485f6488_e81f1340_821.fastq.gz
...
├── barcode06
│   ├── FAW31148_pass_barcode06_11f103b9_993a7465_30.fastq.gz
│   ├── FAW31148_pass_barcode06_11f103b9_993a7465_31.fastq.gz
...
```

The names of subdirectories (`barcode05` and `barcode06`) define the names of samples. All `.fastq.gz` files in each subdirectory are assigned to that sample name.

## Configuration

Run parameters are passed to the workflow by a file called `config.json` that should have these keys:

- `platform`. Either `minion` or `miseq`.
- `samples`. A list containing the sample IDs to analyse.
- `pair`. A list containing one or more of `paired`, `unpaired` and `combined`. This controls whether to analyse paired, unpaired and/or paired and unpaired reads combined. For MinION data `pair` must be `longread`. The concept of paired and unpaired reads does not apply to MinION data but for implementation ease `longread` is used as a placeholder in the workflow wherever `paired`/`unpaired`/`combined` would be used in a MiSeq analysis.
- `order`. A list containing one or both of `primary` and `secondary`. See the IRMA documentation for the distinction between [primary and secondary data](https://wonder.cdc.gov/amd/flu/irma/primary.html) and [residual and secondary assemblies](https://wonder.cdc.gov/amd/flu/irma/secondary_residual.html). If both `primary` and `secondary` are supplied the workflow will do both, and runtime will approximately double.
- `errors`. Either `warn` or `raise`. This controls error handling for some steps in the workflow. `warn` issues warnings if something goes wrong, but will attempt to carry on. `raise` would not carry on.
- `qsr`. A list that can either be empty or contain some combination of `tensqr` and `abayesqr` to
  run a particular quasispecies spectrum reconstruction (QSR) algorithm. See below for more details.

MiSeq example:

```
{
  "platform": "miseq",
  "samples": [
    "YK_2837",
    "YK_2970"
  ],
  "pair": [
    "combined"
  ],
  "order": [
    "primary",
    "secondary"
  ],
  "errors": "warn",
  "qsr": ["tensqr"]
}
```

MinION example:

```
{
  "platform": "minion",
  "samples": [
    "barcode05",
    "barcode06"
  ],
  "pair": [
    "longread"
  ],
  "order": [
    "primary",
  ],
  "errors": "warn",
  "qsr": ["tensqr"]
}
```

(To format JSON nicely do: `python -m json.tool < dirty.json > clean.json`.)

## Preprocessing

Reads are first preprocessed. For MiSeq data adapters are trimmed and quality control reports are generated. For MinION data reads are filtered by a min and max length. See `workflow/rules/preprocess-{minion,miseq}.smk` for details.

```bash
snakemake -s workflow/preprocess.smk -c all
```

`-c all` tells snakemake to use all available cores, scale this back if need be. HTML QC reports for MinION data are saved in `results/qc-raw` and `results/qc-trimmed`.

## Variant calling

Run IRMA and make summary reports of the output:

```bash
snakemake -s workflow/irma.smk -c all
```

### Variant summary reports

Three summary files are generated:

- `results/<order>/xlsx/variants-mcc-by-sample-ordered.xlsx`. In this version,
  each _sample_ has its own sheet, and each sheet contains all influenza
  segments found in that sample.
- `results/<order>/xlsx/variants-mcc-by-segment-ordered.xlsx`. Like above, but
  each _segment_ gets its own sheet.
- `results/<order>/xlsx/variants-mcc-flat-ordered.xlsx`. This contains all
  variants in one flat sheet.

### Sequences

IRMA consensus sequences and amino acid translations are put in
`results/<order>/seq`.

### Splice variants

Splice variants of MP, PA and PB1 are all based on the assumption that IRMA
finds canonical length consensus sequences for these segments (see
[here](splice-variants.md) for more details). 

If IRMA finds a consensus sequence for one of these segments that is not the
expected length, then the behaviour is determined by the config file:

- If `"errors": "warn"` is used, a logfile is produced at the end of the run
(`logs/incorrect-splice-vars-<order>.log`) detailing any segments where the
consensus was not the expected length. For these, it is highly likely that the
variant analysis will be incorrect.
- Instead, if `"errors": "raise"` is used in the config, the workflow stops
immediately if a length mismatch occurs.

(For NS, we know to expect more variability in segment length, and the locations
of exons are flexibly determined.)

## Quasispecies and mixed sample sequence reconstruction

The workflow implements a couple of different programs that try to estimate sequences of segments
that may differ slightly from the reference sequences by either being part of the viral
quasispecies, or from mixed samples. See the [TenSQR](https://github.com/davipatti/TenSQR) and
[aBayesQR](https://github.com/SoYeonA/aBayesQR/) repos and papers for more details on these methods.
Both programs are quite error prone and may fail for unknown reasons on perfectly good input files.
I have tried to put in more catches for the TenSQR errors as it seemed to be [slightly better than
aBayesQR](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293574/).

TenSQR and aBayesQR fail to run when passed BAM files generated by IRMA. So, minimap2, is used to
map all reads (paired and unpaired) to the IRMA combined reference sequence. Therefore, if you want
to run these algorithms then you must have `combined` in the `pair` config file list. 

The inferred sequences get put in `results/qsr/{sample}/{sample}_tensqr.fasta` (or
`{sample}_abayesqr.fasta`), which also contains their estimated frequency in the sample.

If a segment is missing from that file (or that file does not exist) then check:

- `results/qsr/{sample}/{segment}/.tensqr_extractmatrix_log.txt` and
- `results/qsr/{sample}/{segment}/.tensqr_log.txt`

The extract matrix log file comes from the first step in the sequence reconstruction where single
nucleotide variants in BAM files get represented as a binary matrix.


## Bonus: sorted BAM files

Most software to look at reads alignments require sorted bam files, and/or bam
index files. I've written a small workflow for generating these for all bam
files IRMA generates. It requires [samtools](http://www.htslib.org/).
`.sorted.bam` and `.sorted.bam.bai` files are saved in the same directory as the
original `.bam` files. Do:

```bash
snakemake -s workflow/sort-bam.smk -c all
```
