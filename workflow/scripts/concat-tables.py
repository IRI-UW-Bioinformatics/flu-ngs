#!/usr/bin/env python3

import sys
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "concat-tables.py",
        description="Concatenate tables.",
    )
    parser.add_argument("tables", help="TSV files", nargs="+")
    args = parser.parse_args()

    tables = [pd.read_table(path) for path in args.tables]

    pd.concat(tables).to_csv(sys.stdout, sep="\t", index=False)
