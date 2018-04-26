# lacrosse: LAtent ChaRacteristics Of Small-molecule SEnsitivity

A Bayesian nonparametric sparse factor regression for predicting sensitivity of cancer cell lines to a range of drugs. 

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
