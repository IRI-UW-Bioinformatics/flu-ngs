#!/usr/bin/env python3

import argparse
import pandas as pd
from math import inf

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "alter-column-order.py",
        description="Change the order of columns in an excel file.",
    )
    parser.add_argument("--input", help="Input excel file.")
    parser.add_argument("--output", help="Output excel file.")
    parser.add_argument("--order", help="Order of columns.", nargs="+")
    args = parser.parse_args()

    input_xl = pd.ExcelFile(args.input)

    order = {value: index for index, value in enumerate(args.order)}

    with pd.ExcelWriter(args.output) as writer:
        for sheet in input_xl.sheet_names:
            df = input_xl.parse(sheet)
            new_order = sorted(df.columns, key=lambda x: order.get(x, inf))
            df[new_order].to_excel(writer, index=False, sheet_name=sheet)
