#!/usr/bin/env python
"""
Generate a cleaned epitope dataset (HV-V1)

1. Deduplicate full-protein fasta file and convert to dataframe

Since seqs must be unique, we aggregate focal_protein_id into a list
and generate a unique ID (job_name) for each protein sequence.

Many of the proteins are polyproteins- we attempt to filter these out to avoid dealing with 'in-silico cleaving'
and because AF3 has a max token length of ~5120 (https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)

Change `MAX_PROTEIN_LENGTH` to adjust hard threshold.

2. Filter 30mer dataset into non-epitopes and epitopes using the following rules:

- Epitope peptides must have at least a z-score of >= EPITOPE_Z_SCORE_THRESH across at least EPITOPE_MIN_NUM_DONORS_THRESH donors
- Non-epitopes must have a z-score <= NON_EPITOPE_Z_SCORE_THRESH for all donors

3. Take the non-epitope, epitope, and focal protein datasets, and filter them down such that:

- Non-epitope peptides have no overlap with epitope peptides (since, i.e. HV1 contains all peptides,
    reactive or not)
- Focal proteins only contains proteins with >= 1 30-mer from both non-epitope and epitope (as a control)
- Non-epitope only contains peptides which are present in focal proteins after step 2
- Epitope only contains peptides which are present in focal proteins after step 2

Then, annotate non-epitope and epitope with a list column that contains the seq_ids
of the associated focal proteins and the indices into those sequences

Finally, add segmentation boolean masks and localization lists to focal proteins using these annotations.
"""
import polars as pl
import numpy as np
import argparse
from epident.utils import fasta_to_polars, generate_job_name
from pathlib import Path

## Step 1
MAX_PROTEIN_LENGTH = 1500

VALID_AA = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def process_fasta_into_df(fasta_path):
    protein_df = (
        fasta_to_polars(fasta_path)
        .with_columns(
            pl.col("name").str.split("=").list.get(1).alias("raw_protein_id"),
        )
        .select("raw_protein_id", "seq")
    )

    protein_df = protein_df.group_by("seq").agg(
        pl.col("raw_protein_id").alias("raw_protein_ids"),
    )

    protein_df = protein_df.filter(pl.col("seq").str.len_chars() <= MAX_PROTEIN_LENGTH)

    protein_df = protein_df.filter(
        ~pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")
    )

    protein_df = (
        generate_job_name(protein_df, ["seq"], name="job_name")
        .select("job_name", "seq", "raw_protein_ids")
        .sort(by="job_name")
    )
    return protein_df


## Step 2
EPITOPE_Z_SCORE_THRESH = 10
EPITOPE_MIN_NUM_DONORS_THRESH = 2
NON_EPITOPE_Z_SCORE_THRESH = 7

NULL_SPECIES_PLACEHOLDER = 99999999


def process_epitope_data(metadata, z_score):

    donor_colnames = z_score.select(pl.exclude("Sequence name")).columns

    z_score = z_score.with_columns(pl.concat_list(donor_colnames).alias("z_score_list"))

    keep_peptides = metadata

    # add z score lists to peptides
    keep_peptides = keep_peptides.join(
        z_score, left_on="CodeName", right_on="Sequence name"
    )
    keep_peptides = keep_peptides.with_columns(
        pl.col("z_score_list").list.drop_nulls().alias("focal_z_score_list")
    )

    keep_peptides = keep_peptides.with_columns(
        pl.col("focal_z_score_list").list.filter(
            (pl.element().is_not_null()) & (pl.element() >= EPITOPE_Z_SCORE_THRESH)
        )
    )

    keep_peptides = keep_peptides.with_columns(
        pl.col("focal_z_score_list").list.len().alias("focal_z_meet_thresh"),
        pl.col("focal_z_score_list").list.mean().alias("focal_z_mean"),
        pl.col("z_score_list").list.drop_nulls().list.max().alias("all_z_max"),
    )

    non_epitopes = keep_peptides.filter(
        pl.col("all_z_max") < NON_EPITOPE_Z_SCORE_THRESH
    ).with_columns(pl.lit(False).alias("epitope"))

    epitopes = keep_peptides.filter(
        pl.col("focal_z_meet_thresh") >= EPITOPE_MIN_NUM_DONORS_THRESH
    ).with_columns(pl.lit(True).alias("epitope"))

    out_df = pl.concat([non_epitopes, epitopes]).with_columns(
        pl.col("CodeName").alias("raw_peptide_id"),
        pl.col("Peptide").alias("peptide"),
    )

    out_df = (
        generate_job_name(out_df, ["peptide"], name="job_name")
        .select("job_name", "peptide", "raw_peptide_id", "epitope")
        .sort(by="raw_peptide_id")
    )
    return out_df


def convert_null_species_to_int(df):
    df = df.with_columns(
        pl.when(pl.col("SpeciesID").is_not_null())
        .then(pl.col("SpeciesID"))
        .otherwise(pl.lit(NULL_SPECIES_PLACEHOLDER))
    )
    return df


## 3.


def set_30mer_indices_to_true(row_struct):
    """
    Sets specified indices in a boolean mask to True.

    Args:
        row_struct: A dictionary-like object representing a row, with keys
                    "boolmask" and "indices".

    Returns:
        A Polars Series containing the modified boolean mask.
    """
    boolmask = np.array(row_struct["boolmask"])
    indices = row_struct["indices"]
    for idx in indices:
        boolmask[idx : idx + 30] = True
    return boolmask.tolist()


def segment_boolmask_to_localization_list(seg_boolmask):

    seg_boolmask_np = np.array(seg_boolmask)

    # 1D island identification
    in_island = False
    start_idx = None
    islands = []

    for i in range(len(seg_boolmask_np)):

        # start new island
        if not in_island and seg_boolmask_np[i]:
            in_island = True
            start_idx = i

        # end island
        elif in_island and (
            (not seg_boolmask_np[i]) or (i == (len(seg_boolmask_np) - 1))
        ):
            in_island = False
            end_idx = i if not seg_boolmask_np[i] else i + 1
            islands.append([start_idx, end_idx])

    return islands


def annot_segmentation_and_localization(match_df):

    match_df = (
        match_df.with_columns(
            pl.col("seq")
            .str.split("")
            # to get boolmask with equal length to seq
            .list.eval(pl.element() == "NOT_AN_AMINO_ACID")
            .alias("epitope_30mer_segment_boolmask")
        )
        .with_columns(
            pl.col("e_seq_idxs")
            # flattens the list of indices
            .list.eval(pl.element().explode().drop_nulls()).alias("flat_e_seq_idxs"),
        )
        .with_columns(
            pl.struct(
                boolmask="epitope_30mer_segment_boolmask",
                indices="flat_e_seq_idxs",
            )
            .map_elements(set_30mer_indices_to_true, return_dtype=pl.List(pl.Boolean))
            .alias("epitope_30mer_segment_boolmask")
        )
    )

    match_df = match_df.with_columns(
        pl.col("epitope_30mer_segment_boolmask")
        .map_elements(
            segment_boolmask_to_localization_list,
            return_dtype=pl.List(pl.List(pl.Int64)),
        )
        .alias("epitope_localization_list")
    )

    return match_df


def process_epitopes_and_proteins(e_df, ne_df, fp_df):
    match_df_ne = (
        # for those fp sequences that contain a peptide in epitope set
        fp_df.filter(
            pl.col("seq").str.contains_any(
                ne_df.select("peptide").to_series().implode()
            )
        )
        # record all peptides each sequence contains
        .with_columns(
            pl.col("seq")
            .str.extract_many(
                ne_df.select("peptide").to_series().implode(), overlapping=True
            )
            .alias("ne_peptides")
        )
        .explode("ne_peptides")
        .rename({"ne_peptides": "ne_peptide"})
        ## temporarily convert ne_peptide to a 1-length list to allow use of find_many
        .with_columns(pl.concat_list(pl.col("ne_peptide")))
        # and find their indices
        .with_columns(
            pl.col("seq")
            .str.find_many(pl.col("ne_peptide"), overlapping=True)
            .alias("ne_seq_idx")
        )
        ## now convert tmp list back to string
        .with_columns(pl.col("ne_peptide").list.first().alias("ne_peptide"))
        .group_by("job_name", "seq")
        .agg(
            pl.col("ne_peptide").alias("ne_peptides"),
            pl.col("ne_seq_idx").alias("ne_seq_idxs"),
        )
    )

    match_df_e = (
        fp_df.filter(
            pl.col("seq").str.contains_any(e_df.select("peptide").to_series().implode())
        )
        .with_columns(
            pl.col("seq")
            .str.extract_many(
                e_df.select("peptide").to_series().implode(), overlapping=True
            )
            .alias("e_peptides")
        )
        .explode("e_peptides")
        .rename({"e_peptides": "e_peptide"})
        .with_columns(pl.concat_list(pl.col("e_peptide")))
        .with_columns(
            pl.col("seq")
            .str.find_many(pl.col("e_peptide"), overlapping=True)
            .alias("e_seq_idx")
        )
        .with_columns(pl.col("e_peptide").list.first().alias("e_peptide"))
        .group_by("job_name", "seq")
        .agg(
            pl.col("e_peptide").alias("e_peptides"),
            pl.col("e_seq_idx").alias("e_seq_idxs"),
        )
    )

    match_df_all = match_df_e.join(match_df_ne, on=["job_name", "seq"])

    match_df_all = annot_segmentation_and_localization(match_df_all)

    fp_ne_only_df = fp_df.join(match_df_ne.select("job_name"), on="job_name").join(
        match_df_e.select("job_name"), on="job_name", how="anti"
    )

    # filter focal proteins down to only the proteins that contain >1 epitope and >=1 nonepitope
    # this is a control
    fp_df = fp_df.join(
        match_df_all.select(
            "job_name", "epitope_30mer_segment_boolmask", "epitope_localization_list"
        ),
        on="job_name",
    ).sort(by="job_name")

    # annotate epitope and nonepitope df with the job_names and indices into those job_names
    e_annot_filt = (
        match_df_all.select("job_name", "seq", "e_peptides", "e_seq_idxs")
        .explode("e_peptides", "e_seq_idxs")
        .rename(
            {
                "e_peptides": "peptide",
                "job_name": "fp_job_names",
                "e_seq_idxs": "fp_seq_idxs",
            }
        )
        .group_by("peptide")
        .agg(pl.col("fp_job_names"), pl.col("fp_seq_idxs"))
    )

    e_df = e_df.join(e_annot_filt, on="peptide")

    ne_annot_filt = (
        match_df_all.select("job_name", "seq", "ne_peptides", "ne_seq_idxs")
        .explode("ne_peptides", "ne_seq_idxs")
        .rename(
            {
                "ne_peptides": "peptide",
                "job_name": "fp_job_names",
                "ne_seq_idxs": "fp_seq_idxs",
            }
        )
        .group_by("peptide")
        .agg(pl.col("fp_job_names"), pl.col("fp_seq_idxs"))
    )

    ne_df = ne_df.join(ne_annot_filt, on="peptide")

    p_df = pl.concat([e_df, ne_df], how="vertical").sort(by="raw_peptide_id")

    return p_df, fp_df, fp_ne_only_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta_path",
        type=str,
    )
    parser.add_argument("--epitope_data_dir")
    parser.add_argument(
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    output_path = Path(args.output_path)

    # 1.
    protein_df = process_fasta_into_df(args.fasta_path)

    # 2.
    root_data_path = Path(args.epitope_data_dir)

    hv1_metadata = pl.read_csv(
        root_data_path / "PV1_meta_2020-11-23.tsv",
        separator="\t",
    ).filter(pl.col("Category") == "SetCover")

    z_score = pl.read_csv(
        root_data_path / "SHERC_combined_wSB_6-23-21_Z-HDI95.tsv", separator="\t"
    )

    epitope_dat = process_epitope_data(hv1_metadata, z_score)

    # 3.
    e_df = epitope_dat.filter(pl.col("epitope"))
    ne_df = epitope_dat.filter(~pl.col("epitope"))
    ne_df = ne_df.join(e_df.select("peptide"), on="peptide", how="anti")

    epitope_df, full_protein_df, ne_only_protein_df = process_epitopes_and_proteins(
        e_df, ne_df, protein_df
    )

    ne_only_protein_df.write_parquet(output_path / "hv_v1_ne_proteins.parquet")
    epitope_df.write_parquet(output_path / "hv_v1_epitopes.parquet")
    full_protein_df.write_parquet(output_path / "hv_v1_proteins.parquet")
