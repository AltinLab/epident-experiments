#!/usr/bin/env python
"""
Modified from: https://github.com/UberClifford/BepiPred3.0-Predictor/blob/main/bepipred3_CLI.py
"""
### IMPORTS AND STATIC PATHS ###
from bp3 import bepipred3
from pathlib import Path
import argparse
import zipfile
import esm
import torch


class AntigensCached(bepipred3.Antigens):
    """
    Hack the original class so that it can use cached ESM embeddings
    Change the naming scheme of the embeddings file to not include an arbitrary ID.
    Now it is just an accession. Caveat: accessions must now be valid as filenames.
    """

    def get_esm2_represention_on_accs_seqs(self):
        """
        data: list of tuples: [(seq_name, sequence)...]
        per_res_representations:
        """

        # Load ESM-2 model
        if self.run_esm_model_local is not None:
            print(f"Loading ESM2 from: {self.run_esm_model_local}")
            model, alphabet = esm.pretrained.load_model_and_alphabet(
                self.run_esm_model_local
            )

        else:
            model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        batch_converter = alphabet.get_batch_converter()
        model.eval()

        # esm_representations = []
        # preparing batch for ESM2
        upper_case_sequences = [s.upper() for s in self.seqs]
        data = list(zip(self.accs, upper_case_sequences))

        cached_data = []
        new_data = []
        for acc, seq in data:
            if (self.esm_encoding_dir / f"{acc}.pt").is_file():
                cached_data.append((acc, seq, (self.esm_encoding_dir / f"{acc}.pt")))
            else:
                new_data.append((acc, seq))

        print(f"Found {len(cached_data)} cached ESM-2 encodings")

        nr_uncached_seqs = len(new_data)
        batch_generator = self.tuple_generator(new_data)

        enc_id = 0
        completed_new_data = []

        for b in batch_generator:
            batch_labels, batch_strs, batch_tokens = batch_converter(b)
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
            acc_names = batch_labels

            # Extract per-residue representations
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33])
            token_representations = results["representations"][33]

            # encoding sequences
            for i, tokens_len in enumerate(batch_lens):

                esm_representation = token_representations[i, 1 : tokens_len - 1]
                if self.add_seq_len:
                    esm_representation = self.add_seq_len_feature(esm_representation)

                enc_path = self.esm_encoding_dir / f"{acc_names[i]}.pt"
                torch.save(esm_representation, enc_path)

                completed_new_data.append((acc_names[i], batch_strs[i], enc_path))
                enc_id += 1

                print(
                    f"ESM-2 encoded sequence {acc_names[i]} {enc_id}/{nr_uncached_seqs}"
                )

        all_data = cached_data + completed_new_data

        # overwrite these since the order has now changed
        self.accs = [t[0] for t in all_data]
        self.seqs = [t[1] for t in all_data]

        return [t[2] for t in all_data]


### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Make B-cell epitope predictions from fasta file.")
parser.add_argument(
    "-i",
    required=True,
    action="store",
    dest="fasta_file",
    type=Path,
    help="Fasta file contianing antigens",
)
parser.add_argument(
    "-o",
    required=True,
    action="store",
    dest="out_dir",
    type=Path,
    help="Output directory to store B-cell epitope predictions.",
)
parser.add_argument(
    "-pred",
    action="store",
    choices=["mjv_pred", "vt_pred"],
    required=True,
    dest="pred",
    help="Majorty vote ensemble prediction or variable threshold predicition on average ensemble posistive probabilities. ",
)
parser.add_argument(
    "-add_seq_len",
    action="store_true",
    dest="add_seq_len",
    help="Add sequence lengths to esm-encodings. Default is false. On the web server this option is set to true.",
)
parser.add_argument(
    "-esm_dir",
    action="store",
    default=".",
    dest="esm_dir",
    type=Path,
    help="Directory to save esm encodings to. Default is current working directory.",
)
parser.add_argument(
    "-t",
    action="store",
    default=0.1512,
    type=float,
    dest="var_threshold",
    help="Threshold to use, when making predictions on average ensemble positive probability outputs. Default is 0.15.",
)
parser.add_argument(
    "-top",
    action="store",
    default=0.2,
    type=float,
    dest="top_cands",
    help="Top percentage of epitope residues Default is top 20 pct.",
)
parser.add_argument(
    "-rolling_window_size",
    default=9,
    type=int,
    dest="rolling_window_size",
    help="Window size to use for rolling average on B-cell epitope probability scores. Default is 9.",
)
parser.add_argument(
    "-plot_linear_epitope_scores",
    action="store_true",
    dest="plot_linear_epitope_scores",
    help="Use linear B-cell epitope probability scores for plot. Default is false.",
)
parser.add_argument(
    "-z",
    action="store_true",
    dest="zip_results",
    help="Specify option to create zip the bepipred-3.0 results (except the interactive .html figure). Default is false.",
)

args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
var_threshold = args.var_threshold
pred = args.pred
add_seq_len = args.add_seq_len
esm_dir = args.esm_dir
top_cands = args.top_cands
rolling_window_size = args.rolling_window_size
plot_linear_epitope_scores = args.plot_linear_epitope_scores
zip_results = args.zip_results


### FUCNCTIONS ###
def zip_function(result_files, outfile):
    zipf = zipfile.ZipFile(outfile, "w")
    for result_file in result_files:
        zipf.write(result_file, arcname=result_file.name)
    zipf.close()


### MAIN ###

## Load antigen input and create ESM-2 encodings ##

# if you have the esm2 model stored locally, you can this command. To work you need both esm2_t33_650M_UR50D.pt and the esm2_t33_650M_UR50D-contact-regression.pt stored in same directory.
# MyAntigens = bepipred3.Antigens(fasta_file, esm_dir, add_seq_len=add_seq_len, run_esm_model_local=str(WORK_DIR / "models" / "esm2_t33_650M_UR50D.pt") )

MyAntigens = AntigensCached(fasta_file, esm_dir, add_seq_len=add_seq_len)
MyBP3EnsemblePredict = bepipred3.BP3EnsemblePredict(
    MyAntigens, rolling_window_size=rolling_window_size, top_pred_pct=top_cands
)
MyBP3EnsemblePredict.run_bp3_ensemble()

MyBP3EnsemblePredict.create_toppct_files(out_dir)
MyBP3EnsemblePredict.create_csvfile(out_dir)

## B-cell epitope predictions ##
# if pred == "mjv_pred":
#     MyBP3EnsemblePredict.bp3_pred_majority_vote(out_dir)
# elif pred == "vt_pred":
#     MyBP3EnsemblePredict.bp3_pred_variable_threshold(
#         out_dir, var_threshold=var_threshold
#     )

# generate plots (generating graphs for a maximum of 40 proteins)
# MyBP3EnsemblePredict.bp3_generate_plots(
#     out_dir, num_interactive_figs=50, use_rolling_mean=plot_linear_epitope_scores
# )

# zip results
# if zip_results:
#     print("Zipping results")
#     result_files = [f for f in out_dir.glob("*") if f.suffix != ".html"]
#     zip_function(result_files, out_dir / "bepipred3_results.zip")
#     for result_file in result_files:
#         result_file.unlink()
