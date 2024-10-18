#!/usr/bin/env python3

import sys
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "concat-tables.py",
        description="Concatenate TSV data. If all tables are empty, then no file is written.",
    )
    parser.add_argument("tables", help="TSV files", nargs="+")
    args = parser.parse_args()

    tables = []

    for path in args.tables:
        try:
            df = pd.read_table(path)
        except pd.errors.EmptyDataError:
            continue
        else:
            tables.append(df)

    if tables:
        pd.concat(tables).to_csv(sys.stdout, sep="\t", index=False)
