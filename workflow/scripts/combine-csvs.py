#!/usr/bin/env python3

from os.path import basename
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "combine-csvs.py",
        description="""
            Combine multiple CSV files into an excel file with multiple sheets.
        """,
    )
    parser.add_argument("csvs", help="CSV files", nargs="+")
    parser.add_argument("--excel", help="Name of output XLSX file.")
    args = parser.parse_args()

    with pd.ExcelWriter(args.excel) as writer:
        for csv in sorted(args.csvs):
            pd.read_csv(csv).to_excel(writer, sheet_name=basename(csv[:-4]))
