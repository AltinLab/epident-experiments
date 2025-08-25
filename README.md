# Epitope identification

## Overall goals:

[PepSeq](https://www.nature.com/articles/s41596-022-00766-8) is a powerful assay capable of identifying antibody reactivity against a large number of short peptides at a time (~200,000) given a donor blood serum sample and set of peptides. This allows comprehensive diagnostic applications when the peptides come from a group of pathogenic proteins. However, a current limitation of the assay is that only linear epitopes are represented- the short (currently, roughly 30-mer) peptides don't represent the more complex and common conformational epitopes present on the pathogenic proteins. 

As a result, measurements of reactivity are limited- it's possible that a donor is reactive against a pathogenic protein that was tiled into 30-mers for the assay, but that the reactivity couldn't be measured because the donor's antibodies only recognize a conformational epitope from that pathogen.

To overcome this limitation, we can use [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) to scaffold slices of pathogenic proteins into short sequences (roughly 100-mers) still compatible with PepSeq but long enough to fold and recapitulate the epitope region present on the full-length protein.

However, this leaves the problem of creating effective slices of massive libraries of pathogenic proteins that are likely to be conformational epitopes. For example, on a large protein, internally-folded residues are unlikely to contain a conformational epitope region. Naively taking all possible slices of pathogenic proteins is infeasible since the RFdiffusion scaffolding step is highly lossy and computationally expensive- many potential epitopes won't be successfully scaffolded. 

Therefore, we need a model that takes an input protein and outputs a list of high-likelihood conformational epitope targets which can be scaffolded and used in pepseq to measure donor reactivity.

With this ultimate goal in mind, the goals of this epitope identification project are twofold:

### 1. Establish a baseline of epitope identification performance from predicted structure

[BP3](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4497), the currect SOTA model in epitope identification, use ESM-2 LM embeddings to predict where epitope residues are located given a protein sequence. 

While this has proved to be effective, the model's performance was found to increase slightly by adding an RSA (relative solvent accessibility) parameter to the model. The calculation in the BP3  paper comes from [NetSurfP-3.0](https://services.healthtech.dtu.dk/services/NetSurfP-3.0/), which itself uses LM embeddings to compute an RSA estimate.

The RSA parameter was never included in the final, released model due to the slow speed of calculating it. However, the performance boost fromm including the RSA feature, even when approximated without an explicit structure prediction, suggest that structure-derived features are likely to be effective in identifying epitope regions.

The first goal of this project is to re-train a BP3-like architecture using AF3's structural embeddings (some combination of single- and pairwise- embeddings) to attempt to beat BP3's LM-embedding based model on the BP3C50ID test set.

### 2. Aggregate predictions into distinct conformational epitopes

BP3 predicts whether individual residues are a part of an epitope, including information about surrounding residues by using ESM-2 embeddings, which results in a prediction value per each residue in the input protein. This means that the model doesn't identify epitope regions, just epitope residues. This model leaves the problem of aggregation into distinct conformational or linear epitopes unsolved.

The second goal of this project is to aggregate the results of the model from the first stage into distinct conformational epitopes. It may be more effective to train a model explicity for this purpose rather than aggregating a BP3-like model's results, but this will require the creation of a new training and testing dataset with explicit labels for distinct epitopes.