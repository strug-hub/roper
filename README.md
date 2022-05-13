
<!-- README.md is generated from README.Rmd. Please edit that file -->

# roper

<!-- badges: start -->
<!-- badges: end -->

The goal of roper is to accurately detect differentially expressed genes
in the presence of unknown forms of additional variation that remains in
the RNA-Seq read count.

## Installation

You can install the development version of roper from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("strug-hub/roper")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(roper)
set.seed(123)
count <- matrix(rnbinom(n = 5e4, mu = 100, size = 1 / 0.5), ncol = 50)
mod <- model.matrix(~ gl(n = 2, k = 25))
res <- rope(datmat = count, X_model = mod, x_PI_idx = dim(mod)[2])
head(res)
#>          logFC      logLR        adj  adj_logLR      pvals      padj
#> V1 -0.05653202  1.9503623 0.03048107 0.05944913 0.73023249 0.9540478
#> V2  0.03416481  0.6633128 0.01677048 0.01112408 0.88142885 0.9729925
#> V3  0.33041054 70.6842968 0.01985048 1.40311725 0.09389854 0.8661218
#> V4  0.07950771  3.4408076 0.01995857 0.06867361 0.71093261 0.9540478
#> V5  0.06109979  2.4692544 0.01694686 0.04184610 0.77235484 0.9602489
#> V6 -0.06933239  3.6486216 0.01624925 0.05928738 0.73058545 0.9540478
```
