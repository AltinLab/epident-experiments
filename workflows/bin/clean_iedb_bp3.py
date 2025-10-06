#!/usr/bin/env python
"""
Generate the base parquet file from BP3 IEDB paper data
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
        "--raw_data_path",
        type=str,
    )
    parser.add_argument(
        "--discard_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    root_data_path = Path(args.raw_data_path)

    iedb_df = fasta_to_polars(root_data_path / "IEDB_epitopes.fasta").rename(
        {"name": "raw_protein_id"}
    )

    # annotate with the reduced dataset
    iedb_reduced_df = (
        fasta_to_polars(root_data_path / "IEDBSeqsNotSharedAt20ID.fasta")
        .rename({"name": "raw_protein_id"})
        .with_columns(pl.lit(True).alias("reduced"))
    )

    iedb_df = iedb_df.join(
        iedb_reduced_df.select("raw_protein_id", "reduced"),
        on="raw_protein_id",
        how="left",
    ).with_columns(
        pl.col("reduced").fill_null(False).alias("reduced"),
    )

    iedb_df = iedb_df.with_columns(
        pl.col("seq")
        .str.split("")
        .list.eval(pl.element().str.contains(r"[ACDEFGHIKLMNPQRSTVWY]"))
        .alias("epitope_boolmask")
    ).with_columns(pl.col("seq").str.to_uppercase().alias("seq"))

    iedb_df = generate_job_name(
        iedb_df,
        ["seq"],
        name="job_name",
    ).select(
        "job_name",
        "seq",
        "reduced",
        "epitope_boolmask",
        "raw_protein_id",
    )

    bad_seq = iedb_df.filter(
        pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")
    ).with_columns(pl.lit("invalid_characters").alias("reason"))

    if bad_seq.height > 0:
        print(f"Found {bad_seq.height} sequences with invalid characters\n {bad_seq}")

    # filter out bad seq
    iedb_df = iedb_df.filter(~pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"))

    too_long = iedb_df.filter(
        pl.col("seq").str.len_chars() > MAX_PROTEIN_LENGTH
    ).with_columns(pl.lit("too_long").alias("reason"))

    if too_long.height > 0:
        print(
            f"Found {too_long.height} sequences longer than {MAX_PROTEIN_LENGTH} characters\n {too_long}"
        )

    iedb_df = iedb_df.filter(pl.col("seq").str.len_chars() <= MAX_PROTEIN_LENGTH)

    # sometimes, a sequence appears multiple times
    # we deal with this by checking if the different epitopes overlap
    # if they do, we throw out the sequence, since the epitope region is ambiguous
    # if they don't, we merge them, so the sequence has multiple (possibly contiguous)
    # epitope regions
    multi_row_seqs = (
        # just for determinism
        iedb_df.sort(by="seq")
        .group_by("seq")
        .agg(
            pl.col("epitope_boolmask").alias("epitope_boolmask_list"),
        )
        # these potentially have overlaps
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
        # if so, we will filter it out later
        # since it is ambiguous where the epitope starts and stops
        .with_columns(
            pl.col("epitope_boolmask_list_T")
            .list.eval(pl.element().list.eval(pl.element().sum() > 1).flatten())
            .alias("overlap_per_residue_boolmask"),
            # if not, we will merge all the rows that share a sequence, performing a bitwise AND on their boolmasks
            pl.col("epitope_boolmask_list_T")
            .list.eval(pl.element().list.eval(pl.element().sum() > 0).flatten())
            .alias("merged_epitope_boolmask_list"),
        )
        .with_columns(
            pl.col("overlap_per_residue_boolmask").list.any().alias("overlap")
        )
    )

    epitopes_with_overlap = (
        (
            (
                multi_row_seqs.filter(pl.col("overlap"))
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
            .join(iedb_df, on=["seq", "epitope_boolmask"], how="inner")
        )
        .select(
            [
                "job_name",
                "seq",
                "reduced",
                "epitope_boolmask",
                "raw_protein_id",
            ]
        )
        .with_columns(pl.lit("overlap").alias("reason"))
    )

    epitopes_no_overlap = (
        multi_row_seqs.filter(~pl.col("overlap"))
        .rename(
            {
                "merged_epitope_boolmask_list": "epitope_boolmask",
            }
        )
        .select(
            "seq",
            "epitope_boolmask",
        )
    )

    epitopes_no_overlap = (
        epitopes_no_overlap.join(
            iedb_df.select("job_name", "seq", "raw_protein_id", "reduced"),
            on="seq",
            how="inner",
        )
        .group_by("job_name", "seq", "epitope_boolmask")
        .agg(
            pl.col("reduced").first(),
            pl.col("raw_protein_id").alias("raw_protein_ids"),
        )
    ).select(
        [
            "job_name",
            "seq",
            "reduced",
            "epitope_boolmask",
            "raw_protein_ids",
        ]
    )

    if epitopes_with_overlap.height > 0:
        print(
            f"Found {epitopes_with_overlap.height} sequences with ambiguous epitope regions\n {epitopes_with_overlap}"
        )

    if epitopes_no_overlap.height > 0:
        print(
            f"Found {epitopes_no_overlap.explode(['raw_protein_ids']).height} sequences with multiple epitope regions, combining them\n {epitopes_no_overlap}"
        )

    iedb_df = iedb_df.join(multi_row_seqs.select("seq").unique(), on="seq", how="anti")

    iedb_df = iedb_df.group_by(pl.exclude("raw_protein_id")).agg(
        pl.col("raw_protein_id").alias("raw_protein_ids")
    )

    iedb_df = pl.concat(
        [
            iedb_df,
            epitopes_no_overlap,
        ]
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

    iedb_df.write_parquet(
        args.output_path,
    )
