# HiSSE model for the geographical spread of SARS-CoV-2 

This repository contains code for fitting hidden state multistate speciation and extionsion (HiSSE) model using the SARS-CoV-2 genomes associated with the country of which they originate.


## Reading the source files

If you are a statistician and you are only interested in reading the code for the actual model, please see the six `*.rb` files for the HiSSE model. They are RevBayes model files (similar to Stan, JAGS files) which contains the formulas for the prior and likelihood. For the MuSSE model, please see `full-tree-est.R` (MLE), `full-tree-boot.R` (bootstrap) and `full-tree-bootanalys.R` (computing CI out of the bootstarp results). Other files are mostly plotting scripts. `preprocess_ST6NONULTRA.R` contains all pre-processing code from the raw Next Strain data.


## Data Preprocessing

The following pre-processing steps are done to the phylogeny:

1. Rescale edge length so that a unit length represents 30 days
2. Resolve multichotomies with `ape::multi2di`
3. `ape::multi2di` produces dichotomies with some length-zero branches. We replace zero length with 1 hour.
4. Use `ape::collapse.singles` to resolve internal nodes with a single descendant.


## Numerical results

The results of running the scripts are zipped in `numerical-results.tar.xz` to
stop Github complaining about big file size. Most importantly, the zip file
contains `model.log`, CSV file that contains the MCMC chain, along with some
plots and other files that RevBayes outputs, some of them are needed by the
plotting scripts.


# How to run the experiment

## Requirements

1. [RevBayes >= 1.0.12](https://revbayes.github.io)
2. [GNU R](https://www.r-project.org) or [pqR](http://www.pqr-project.org)


## Data

The raw data from Next Strain is in `DATA_MUSSE2`. There is a pre-processing script `preprocess_ST6NONULTRA.R` which the user should run before running the MCMC estimation. This script transforms the phylogeny by re-scale the tree edges (for numerical stability), collapsing singleton internal nodes, etc. To pre-process the tree, do 

```
mkdir ALL_CLEANED_DAT_ST6NOULTRA
R --vanilla < preprocess_ST6NONULTRA.R
```

Now the new pre-processed, cleaned, oven-ready data is stored in `ALL_CLEANED_DAT_ST6NOULTRA`. In this folder is the version of the data that the MCMC estimation script will read.


## Running the Bayesian MCMC estimation

There are six models, which correspond to six `.rb` files. For example, to run the `ALLEQ` model, please use

```
rb-mpi covid19-genome-hisse-ALLEQ.rb
```

After this finished or interrupted by `Ctrl-C`, the main MCMC result will be stored at `RES_HISSE_ST6NOULTRA_ALLEQ_1/model.log`, which is a CSV file.

## Plotting

To plot the model, the user needs to run `postproc_*.rb`, `rbana_ST6.R` and `plot_simmap4.R` in the respective order. The plots are saved in `RES_HISSE_ST6NOULTRA_*_1/*.pdf`. Instead of `plot_simmap4.R`, you may also try other `plot_simmap*.R`; they produces similar tree plots but with different colouring.
