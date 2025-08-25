#!/usr/bin/env python
"""
Generate the base parquet file from BP3 paper data
"""
import polars as pl
from pathlib import Path
import argparse
from epident.utils import generate_job_name, fasta_to_polars
import numpy as np

MAX_PROTEIN_LENGTH = 1500


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--train_fasta",
        type=str,
    )
    parser.add_argument("--test_fasta", type=str)
    parser.add_argument(
        "--discard_path",
        type=str,
    )
    parser.add_argument(
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    bp_test_df = (
        fasta_to_polars(Path(args.train_fasta))
        .rename({"name": "raw_protein_id"})
        .with_columns(
            pl.lit(True).alias("test"),
        )
    )

    bp_train_df = (
        fasta_to_polars(Path(args.test_fasta))
        .rename({"name": "raw_protein_id"})
        .with_columns(
            pl.lit(False).alias("test"),
        )
    )

    bp_df = pl.concat([bp_test_df, bp_train_df])

    bp_df = bp_df.with_columns(
        pl.col("seq")
        .str.split("")
        .list.eval(pl.element().str.contains(r"[ACDEFGHIKLMNPQRSTVWY]"))
        .alias("epitope_boolmask")
    ).with_columns(pl.col("seq").str.to_uppercase().alias("seq"))

    bp_df = generate_job_name(
        bp_df,
        ["seq"],
        name="job_name",
    ).select(
        "job_name",
        "seq",
        "test",
        "epitope_boolmask",
        "raw_protein_id",
    )

    bad_seq = bp_df.filter(
        pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")
    ).with_columns(pl.lit("invalid_characters").alias("reason"))

    if bad_seq.height > 0:
        print(f"Found {bad_seq.height} sequences with invalid characters\n {bad_seq}")

    # filter out bad seq
    bp_df = bp_df.filter(~pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"))

    too_long = bp_df.filter(
        pl.col("seq").str.len_chars() > MAX_PROTEIN_LENGTH
    ).with_columns(pl.lit("too_long").alias("reason"))

    if too_long.height > 0:
        print(
            f"Found {too_long.height} sequences longer than {MAX_PROTEIN_LENGTH} characters\n {too_long}"
        )

    bp_df = bp_df.filter(pl.col("seq").str.len_chars() <= MAX_PROTEIN_LENGTH)

    epitope_overlaps = (
        # just for determinism
        bp_df.sort(by="seq")
        .group_by("seq")
        .agg(
            pl.col("epitope_boolmask").alias("epitope_boolmask_list"),
        )
        .filter(pl.col("epitope_boolmask_list").list.len() > 1)
        .with_columns(
            pl.col("epitope_boolmask_list")
            .map_elements(
                lambda x: (np.stack(x.to_numpy()).T).tolist(),
                return_dtype=pl.List(pl.List(pl.Boolean)),
            )
            .alias("epitope_boolmask_list_T")
        )
        # is a residue ever called an epitope in more than one sequence?
        .with_columns(
            pl.col("epitope_boolmask_list_T")
            .list.eval(pl.element().list.eval(pl.element().sum() > 1).flatten())
            .alias("overlap_per_residue_boolmask")
        )
        .with_columns(
            pl.col("overlap_per_residue_boolmask").list.any().alias("overlap")
        )
    )

    epitopes_with_overlap = (
        (
            (
                epitope_overlaps.filter(pl.col("overlap"))
                .explode(["epitope_boolmask_list"])
                .rename(
                    {
                        "epitope_boolmask_list": "epitope_boolmask",
                    }
                )
            )
            .select(
                "seq",
                "epitope_boolmask",
            )
            .join(bp_df, on=["seq", "epitope_boolmask"], how="inner")
        )
        .select(
            [
                "job_name",
                "seq",
                "test",
                "epitope_boolmask",
                "raw_protein_id",
            ]
        )
        .with_columns(pl.lit("overlap").alias("reason"))
    )

    bp_df = bp_df.join(
        epitopes_with_overlap.select("job_name"), on=["job_name"], how="anti"
    )

    discard = pl.concat(
        [
            bad_seq,
            too_long,
            epitopes_with_overlap,
        ]
    )

    discard.write_parquet(
        args.discard_path,
    )

    bp_df.write_parquet(
        args.output_path,
    )
