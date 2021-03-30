# HiSSE model for the Geographical Spread of SARS-CoV-2

This repository contains code for fitting hidden state multistate speciation and extionsion (HiSSE) model using the SARS-CoV-2 genomes associated with the country of which they originate.


## Reading the source files

If you are a statistician and you are only interested in reading the code for the actual model, please see the six `*.rb` files for the HiSSE model. They are RevBayes model files (similar to Stan, JAGS files) which contains the formulas for the prior and likelihood. The `mpp-*.ksh` scripts are used for submitting the RebBayes programs to a Slurm job queue for high-performance computation; you don't need them if you don't use Slurm workload manager. `run_me_to_plot_MCMC.ksh` is a simple ksh script which calls `plotMC.R` with the correct file names, so you may want to read the latter to see how graphs are plotted. `plotMC.R` also sources some other .R files which does the actual graph plotting and computes the BIC of each models. For the MuSSE model, which does not converge, please see `musse-mle-withmu.R`, which fits the MuSSE model with conjugate gradient optimisation. 


## Data Preprocessing

The following pre-processing steps are done to the phylogeny:

1. Rescale edge length so that a unit length represents 30 days
2. Resolve multichotomies with `ape::multi2di`
3. `ape::multi2di` produces dichotomies with some length-zero branches. We replace zero length with 1 hour.
4. Use `ape::collapse.singles` to resolve internal nodes with a single descendant.


## Numerical results

The files are currently too large to upload here and it will be uploaded soon,
probably in a form of splitted zip files.


# How to run the experiment

## Requirements

1. [RevBayes == 1.0.12](https://revbayes.github.io)
2. [GNU R](https://www.r-project.org)
3. Some CRAN packages (see source code)
4. ksh93, which is available in almost all Linux distros package managers.

## Data

The raw data from Next Strain is in `RAW_PHYLOGENY`. There is a pre-processing script `preprocess_ST6NONULTRA.R` which the user should run before running the MCMC estimation. This script transforms the phylogeny by re-scale the tree edges (for numerical stability), collapsing singleton internal nodes, etc. To pre-process the tree, do 

```
mkdir ALL_CLEANED_DAT_ST6NOULTRA
R --vanilla < preprocess_ST6NONULTRA.R
```

Now the new pre-processed, cleaned, oven-ready data is stored in `ALL_CLEANED_DAT_ST6NOULTRA`. In this folder is the version of the data that the MCMC estimation script will read.


## Running the Bayesian MCMC estimation

There are six models, which correspond to six `.rb` files. The `mpp-*-withmu.ksh` should illustrate directly how to run these scripts.
As an example, to run the `ALLEQ` model, use the following bash command:

```
    echo \
        'chnname="01"; '   \
        'outdir="OUTPUT_DIR"; ' \
        'source("covid19-genome-hisse-alleq-withmu.rb")'| \
        rb
```

Note that when we run the experiment, we assume that RevBayes was compiled with OpenMPI turnt *off*. Here we use concurrency to run multiple parallel chains rather than using RevBayes' built-in OpenMPI feature. After this has finished or interrupted by `Ctrl-C`, the main MCMC result will be stored at `OUTPUT_DIR/model01.log`, which is a CSV file despite of its file name.


## Plotting

To plot the MCMC chain's trace plot and the marginal posterior, use the `run_me_to_plot_MCMC.ksh`. The script plots the graphs by simply
running a for-loop which evokes `./plotMC.R` with the correct directory name. If you want to plot the estimated HiSSE states, you should
run `./make_charmap.ksh` first (see the file for more detail).

