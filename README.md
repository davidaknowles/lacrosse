# lacrosse: LAtent ChaRacteristics Of Small-molecule SEnsitivity

A Bayesian nonparametric sparse factor regression for predicting sensitivity of cancer cell lines to a range of drugs. 

Drug screening studies assay the sensitivity of a range of cancer cell lines across an array of anti-cancer therapeutics. Alongside these sensitivity measurements high dimensional molecular characterizations of the cell lines are available, including gene expression, copy number variation and genomic mutations. `lacrosse` is a sparse multitask regression model which learns discriminative latent characteristics that predict drug sensitivity and are associated with specific molecular features. We use Bayesian nonparametrics (specifically an Indian buffet process prior) to automatically infer the appropriate number of these latent characteristics. The resulting analysis couples high predictive performance with interpretability since each latent characteristic involves a typically small set of drugs, cell lines and genomic features. `lacrosse` uncovers a number of drug-gene sensitivity associations missed by single gene analyses.

## Installation

The `lacrosse` R package itself can be installed using e.g.
```
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("davidaknowles/lacrosse/lacrosse")
```
which will automatically install the required dependencies. 

To run `lacrosse` on the CTRPv2 dataset you will need the `PharmacoGx` R package [https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv723], installed as
```
source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
```

You'll also need the `tidyverse`, `magrittr` and `doMC` (or `foreach`) packages installed. 

## Usage 

`code/run_lacrosse.R` gives a working example of running `lacrosse` on the CTRPv2 data. 
