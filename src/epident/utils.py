from Bio import SeqIO
import polars as pl
import hashlib


def hash_sequence(seq: str, hash_type: str = "md5") -> str:
    """
    Hash a TCR sequence using the specified hash function.

    Args:
        tcr_seq (str): The TCR sequence string.
        hash_type (str): The hash function to use ('md5', 'sha1', 'sha256', etc.)

    Returns:
        str: The hexadecimal digest of the hashed sequence.
    """
    # Select the hash function
    if hash_type.lower() == "md5":
        h = hashlib.md5()
    elif hash_type.lower() == "sha1":
        h = hashlib.sha1()
    elif hash_type.lower() == "sha256":
        h = hashlib.sha256()
    else:
        raise ValueError("Unsupported hash type")

    # Encode the sequence and compute the hash
    h.update(seq.encode("utf-8"))
    return h.hexdigest()


def fasta_to_polars(fasta_path: str, desc_as_name: bool = False) -> pl.DataFrame:
    """
    Read a FASTA file and convert it into a Polars DataFrame
    with columns ["name", "sequence"].

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    pl.DataFrame
        - "name": the sequence ID
        - "sequence": the full sequence string
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if desc_as_name:
        names = [rec.description for rec in records]
    else:
        names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]

    df = pl.DataFrame({"name": names, "seq": seqs})
    return df


def generate_job_name(df, cols, name="job_name"):
    df = df.with_columns(
        pl.concat_str(
            pl.concat_str(
                [
                    *[pl.col(colname) for colname in cols],
                ],
                ignore_nulls=True,
            )
            .map_elements(lambda x: hash_sequence(x, "md5"), return_dtype=pl.String)
            .alias(name),
        )
    )
    return df
