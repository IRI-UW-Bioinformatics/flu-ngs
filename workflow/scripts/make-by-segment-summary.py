#!/usr/bin/env python3

import argparse
import pandas as pd


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "make-by-segment-summary.py",
        description="""Make a summary excel workbook, where each sheet contains data from
                    different samples, but the same segment. At the same time write a
                    workbook that has all data on a single (flat) sheet.""",
    )
    parser.add_argument("--in-excel", help="By sample excel sheet.")
    parser.add_argument("--out-segment", help="Name of output by-segment excel sheet.")
    parser.add_argument("--out-flat", help="Name of single sheet output.")
    args = parser.parse_args()

    xl = pd.ExcelFile(args.in_excel)
    dfs = []
    for sheet in xl.sheet_names:
        df = xl.parse(sheet)
        df["Sample"] = sheet
        dfs.append(df)

    df = pd.concat(dfs)

    with pd.ExcelWriter(args.out_segment) as writer:
        for segment, sub in df.groupby("Segment"):
            sub.to_excel(writer, sheet_name=segment, index=False)

    df.to_excel(args.out_flat, index=False)
