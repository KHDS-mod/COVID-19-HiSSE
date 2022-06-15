# HiSSE model for the Geographical Spread of SARS-CoV-2

This repository contains code for fitting hidden state multistate speciation and extionsion (HiSSE) model using the SARS-CoV-2 genomes associated with the country of which they originate.


## Reading the source files

If you are a statistician and you are only interested in reading the code for the actual model, please see the six `*.rb` files for the HiSSE model. They are RevBayes model files (similar to Stan, JAGS files) which contains the formulas for the prior and likelihood. The `mpp-*.ksh` scripts are used for submitting the RevBayes programs to a Slurm job queue for high-performance computation; you don't need them if you don't use Slurm workload manager. `plotMC.R` can be sourced into R to produce graphs about the MCMC. The script also sources some other .R files which does the actual graph plotting and computes the BIC of each models. For the MuSSE model, which does not converge, please see `musse-mle-withmu.R`, which fits the MuSSE model with conjugate gradient optimisation. 


## Data Preprocessing

The following pre-processing steps are done to the phylogeny:

1. Rescale edge length so that a unit length represents 30 days
2. Resolve multichotomies with `ape::multi2di`
3. `ape::multi2di` produces dichotomies with some length-zero branches. We replace zero length with 1 hour.
4. Use `ape::collapse.singles` to resolve internal nodes with a single descendant.

## Numerical results (directly loadable into R, without stochastic character maps)

The easiest way to get the MCMC chain loaded into R is to download [this ".rds" file](https://github.com/KHDS-mod/COVID-19-HiSSE/blob/4a3a9fb9cef5998f333c1014d0237ff5f3784d50/numerical_result_monad.rds.lrz), which contains all the chains of all the models, including both the burn-in-trimmed version and the raw version, in a single R object. This R object also contains all the steps and code
that was used to read in, trim and postprocess the MCMC chains. See [here](https://github.com/KHDS-mod/COVID-19-HiSSE/tree/rmonadRDS) for how to decompress a `.lrz` file format. You will need the R package [rmonad](https://cran.r-project.org/web/packages/rmonad/) (0.7.0 was used) to manipulate the object.

Use the following to commands in R to get the chains of, for example, the full model with informative prior:

```
library(rmonad)
library(magrittr)

m = readRDS('numerical_result_monad.rds')
chns = view(m, 'MCMCchains:full-prior3') %>% esc
```

The burn-in-trimmed counterpart can be obtained with

```
chns_trimmed = view(m, 'MCMCchains_trimmed:full-withmu') %>% esc
```

You can get all tags and explore the [rmonad](https://cran.r-project.org/web/packages/rmonad/) object by using this command:

```
m %>% get_tag %>% unlist %>% unique
```


## Numerical results (from raw file, with stochastic character maps)

You can download the raw files of the MCMC chains [here](http://urn.kb.se/resolve?urn=urn%3Anbn%3Ase%3Aliu%3Adiva-185867). Note that this file is about 48GB, much larger than the previously mentioned rmonad object because it contains all the inferred stochastic character maps.

After extracting the compressed file, this project folder should look like the following:

<pre>
COVID-19-HiSSE
<b>├── chain_data
    ├── RES_HISSE_ST6NOULTRA_ALLEQ_WITHMU_2
    │   ├── model01.log
    │   ├── model02.log
    │   ├── model03.log
        ...
    ├── RES_HISSE_ST6NOULTRA_FULL_WITHMU_2
    │   ├── model01.log
    │   ├── model02.log
    │   ├── model03.log
        ...
</b>├── ALL_CLEANED_DAT_ST6NOULTRA
├── ALL_CLEANED_DAT_ULTRAMETRIC
├── AirPax
├── RAW_PHYLOGENY
...
...
...
</pre>

Some models has a "PART1"-suffixed directory because their MCMC chains were stopped and re-started because the computation cluster that we used to run the models on has a hard-limit maximal duration for each process. These "PART1" chains should be concatenated with the chains contained in their counterpart folder that does not has the "PART1" suffix. Running the `plotMC.R` script does this concatenation as well as reproducing the [rmonad](https://cran.r-project.org/web/packages/rmonad/) object mentioned above, which is stored in a global variable called `m`. Before running `plotMC.R` you must make an empty folder called `plots` in the project's root folder.

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

To plot the MCMC chain's trace plot and the marginal posterior, use the `plotMC.R` script. If you want to plot the estimated HiSSE states (the painted tree graph), you should run `./make_charmap.ksh` first (see the file's content for more detail, and you will need to change the model name in the script to get the character map that you want). The script for producing the painted tree graph is in `plot_simmap*.R`.

