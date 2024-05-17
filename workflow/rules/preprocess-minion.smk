"""
Preprocessing sequence data from the MinION platform.
"""

from pathlib import Path


rule combine_filter_reads:
    """Combine multiple MinION reads in a single directory and filter them based on min and max lengths."""
    input:
        lambda wild: Path(f"raw/{wild.sample}").glob("*.fastq.gz"),
    output:
        "processed_reads/{sample}/{sample}_filtered.fastq.gz",
    shell:
        "cat {input} | gunzip -c | chopper -q 10 --minlength 600 --maxlength 2500 | gzip > {output}"
