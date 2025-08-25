#!/usr/bin/env python
import polars as pl
from epident.sasa import PatchedSASAAnalysis
from mdaf3.FeatureExtraction import split_apply_combine
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path
import argparse


def sasa(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    analysis = PatchedSASAAnalysis(u)
    analysis.run()
    row["RSA"] = analysis.results.relative_residue_area[0].tolist()
    row["SA"] = analysis.results.residue_area[0].tolist()

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

    pq = pl.read_parquet(args.input_parquet)
    pq = split_apply_combine(pq, sasa, Path(args.inference_path), chunksize=15)

    pq.write_parquet(args.output_path)
