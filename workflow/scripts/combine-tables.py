#!/usr/bin/env python3

from os.path import basename
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "combine-tables.py",
        description="""
            Combine multiple flat table files into an excel file with multiple sheets.
            Tables could be comma-separated, tab-separated, or anything understood by the
            `pandas.read_table` `sep=` argument.
        """,
    )
    parser.add_argument("tables", help="Table files", nargs="+")
    parser.add_argument("--excel", help="Name of output XLSX file.")
    parser.add_argument("--sep", help="Separator used in the input files.", default="\t")
    parser.add_argument(
        "--comment",
        help="Indicates a comment in the inputs.",
        default=None,
        required=False,
    )
    args = parser.parse_args()

    with pd.ExcelWriter(args.excel) as writer:
        for table in sorted(args.tables):
            name = basename(table)
            name = name[: name.rindex(".")]

            try:
                df = pd.read_table(
                    table, sep=args.sep, comment=args.comment, engine="python"
                )
            except pd.errors.EmptyDataError:
                df = pd.DataFrame({"No data": []})

            df.set_index(df.columns[0]).to_excel(writer, sheet_name=name)
