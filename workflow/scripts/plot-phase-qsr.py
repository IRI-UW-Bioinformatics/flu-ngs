from pathlib import Path
import argparse

import pandas as pd


def add_snp_column(df: pd.DataFrame) -> pd.DataFrame:
    df["is_snp"] = (df["Consensus_Allele"].str.len() == 1) & (
        df["Minority_Allele"].str.len() == 1
    )
    return df


def load_irma_tables(paths: list[str]) -> pd.DataFrame:
    tables = [pd.read_table(table) for table in paths]
    tables = [table for table in tables if not table.empty]
    if tables:
        return pd.concat(tables).pipe(add_snp_column).query("is_snp")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("plot-phase-qsr.py")
    parser.add_argument("--irma-tables", nargs="+")
    args = parser.parse_args()

    df_irma = load_irma_tables(args.irma_tables)

