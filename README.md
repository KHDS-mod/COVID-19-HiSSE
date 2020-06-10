# HiSSE model for the geographical spread of SARS-CoV-2 

This repository contains code for fitting hidden state multistate speciation and extionsion (HiSSE) model using the SARS-CoV-2 genomes associated with the country of which they originate.


## Reading the source files

If you are a statistician and you are only interested in reading the code for the actual model, please see the six `*.rb` files. They are RevBayes model files (similar to Stan, JAGS files) which contains the formulas for the prior and likelihood. Other files are mostly pre-processing and plotting scripts. `preprocess_ST6NONULTRA.R` contains all pre-processing code from the raw Next Strain data.


## Requirements

1. [RevBayes >= 1.0.12](https://revbayes.github.io)
2. [GNU R](https://www.r-project.org) or [pqR](http://www.pqr-project.org)


## Data

The raw data from Next Strain is in `DATA_MUSSE2`. There is a pre-processing script `preprocess_ST6NONULTRA.R` which the user should run before running the MCMC estimation. This script transforms the phylogeny by re-scale the tree edges (for numerical stability), collapsing singleton internal nodes, etc. To pre-process the tree, do 

```
mkdir ALL_CLEANED_DAT_ST6NOULTRA
R --vanilla < preprocess_ST6NONULTRA.R
```

Now the new pre-processed, cleaned, oven-ready data is stored in `ALL_CLEANED_DAT_ST6NOULTR`. In this folder is the version of the data that the MCMC estimation script will read.


## Running the Bayesian MCMC estimation

There are six models, which corresponds to six `.rb` files. For example, to run the `ALLEQ`, please use

```
rb-mpi covid19-genome-hisse-ALLEQ.rb
```

After this finished or the user use `Ctrl-C` to interrupt, the main MCMC result will be stored at `RES_HISSE_ST6NOULTRA_ALLEQ_1/model.log`, which is a CSV file.

## Plotting

To plot the model, the user needs to run `postproc_*.rb`, `rbana_ST6.rb` and `plot_simmap4.R` in the respective order. The plots are saved in `RES_HISSE_ST6NOULTRA_*_1/*.pdf`.
