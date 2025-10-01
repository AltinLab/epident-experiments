#!/usr/bin/env python
from mdaf3.FeatureExtraction import split_apply_combine
from pathlib import Path
import polars as pl
import argparse

from mdaf3.AF3OutputParser import AF3Output


def extract_residue_pLDDT(row, inf_path):
    af3 = AF3Output(inf_path / row["job_name"])
    row["pLDDT"] = [
        float(atm_pLDDT.mean())
        for atm_pLDDT in af3.get_mda_universe().residues.tempfactors
    ]
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument(
        "--inference_path",
        type=str,
    )
    parser.add_argument(
        "--output_path",
        type=str,
    )

    args = parser.parse_args()

    fp = pl.read_parquet(args.input_parquet)

    fp = split_apply_combine(
        fp, extract_residue_pLDDT, Path(args.inference_path), chunksize=15
    )

    fp.write_parquet(args.output_path)
