# Analyses and data from Castledine_ et al. (2018) JEB

Data, analyses and figures from: _Castledine et al. (2018) A shared coevolutionary history does not alter the outcome of coalescence in experimental populations of Pseudomonas fluorescens. Journal of Evolutionary Biology_ 

DOI of paper:



DOI of analyses and dataset:

[![DOI](https://zenodo.org/badge/153314427.svg)](https://zenodo.org/badge/latestdoi/153314427)

### Outline

This repository contains the final datasets, R scripts containing the analyses, and figures of the above-mentioned paper. It can recreate the analyses and figures in the main text (Figure 2 tand Figure 3). Figure 1 is present as a `pdf`.

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Castledine_et_al_2018_JEB/issues) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in `scripts/`.
- `scripts/morph_level_analysis.R` contains the analysis at the morphotype level. This script tests whether genotypes benefit from the presence of another coevolved genotyope. It also looks at whether genotypes benefit from the presence of another genotype of the same strain (i.e. LacZ or wild-type). Recreates Figure 2.
- `scripts/community_level_analysis.R` contains the analysis of trials where both morphs (smooth and wrinkly spreader) were present in both wild-type and LacZ inoculates. This script tests whether coevolution changes the outcome of community coalescence and recreates Figure 3.

### Overview of data files

- column meanings for `data/all_trials_morphs_counts.csv`
    - `comb`: unique code for each competition trial
    - `WT_div`: number of wild-type morphs in the trial
    - `LZ_div`: number of LacZ morphs in the trial
    - `type`: whether the counts are for the wild-type or LacZ strain
    - `morph`: the morph of the specific count
    - `community`: the community of the morph
    - `coevolved`: whether or not the morph had a coevolved morph present in the trial
    - `morph_comb`: the morph combination in the trial for each strain
    - `WT_coev`: whether the wild-type community shared a coevolutionary history
    - `LZ_coev`: whether the LacZ community shared a coevolutionary history
    - `treatment`: the number of LacZ morphs:wild-type morphs and their coevolutionary history
    - `treatment2`: the LacZ morph ID:wild-type morph ID and their coevolutionary history
    - `T0_count_ml`: number of cells inoculated per ml
    - `T1_count`: observed number of CFUs on the plate
    - `T1_count_ml`: final count transformed to count per ml

- column meanings for `data/all_combs.csv`
    - `Combination`: unique code for each competition trial
    - `WT`: input morphs of the wild-type populations and the community they came from
    - `LacZ`: input morphs of the LacZ populations and the community they came from

__All analyses are done in R version 3.5.1, on macOS High Sierra 10.13.6. I am unsure whether some older versions of R will support all of the packages.__
