
<!-- README.md is generated from README.Rmd. Please edit that file -->

# roper

<!-- badges: start -->
<!-- badges: end -->

<<<<<<< HEAD
The goal of roper is to accurately detect differential expressed genes
in the presence of unknown forms of additional variation that remains in
the RNA-Seq read count.
=======
The goal of roper is to accurately detects differential expressed genes
in the presence of unknown forms of additional variation that remains in
the RNA-Seq read count
>>>>>>> 76aeb48eb43942fee25835d7b67ed34037546724

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
count <- matrix(rnbinom(n = 5e4, mu = 100, size = 1 / 0.5), ncol = 50)
mod <- model.matrix(~ gl(n = 2, k = 25))
res <- rope(datmat = count, X_model = mod, x_PI_idx = dim(mod)[2])
```
