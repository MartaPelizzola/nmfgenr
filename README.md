
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nmfgenr

<!-- badges: start -->

<!-- badges: end -->

The goal of nmfgenr is to provide an easy to use tool to perform
traditional and convex non-negative matrix factorization with a wide
range of underlying distributions. This package provides a fast
implementation of updates derived via MM algorithms for traditional and
convex NMF under the following distributions: \* Tweedie (with its
special cases Normal and Poisson) \* Negative Binomial \* Laplace.

## Installation

You can install the development version of nmfgenr as follows:

``` r
#devtools::install_github("MartaPelizzola/nmfgenr",force = TRUE, build_vignettes = TRUE)
```

## Basic example

This is a basic example for the traditional Tweedie NMF with power = 0
(corresponding to the Normal distribution):

``` r
## Load the package 
library(nmfgenr)
## Simulate a small data set to run examples:
data <- matrix(rnorm(100, mean = 50, sd = 10), nrow = 10)
## Run traditional Tweedie NMF:
#res <- nmfgen(data, rank = 3, "Tweedie", "traditional", pwr = 0)
```

The output can then be accessed as follows:

``` r
## Access weight matrix
#head(res$W)
## Access features matrix
#head(res$H)
```

Similarly, all other methods can be run in this way by specifying the
distribution and the method parameters.

For convex NMF the results can be accessed as follows:

``` r
## Access weight matrix
#res$W1
## Access features matrix
#res$W2
```

## Example with custom initial values

We also provide the option to set custom initial values. This can be
done by adding input values for wmat and hmat (if traditional NMF is
used) or for w1mat and h1mat (if convex NMF is used). For example, for
traditional Tweedie NMF with power = 0, it can be done as follows:

``` r
## Create an initial value for the weight matrix:
winit <- matrix(runif(nrow(data)*3), nrow = nrow(data))
## Create an initial value for the feature matrix:
hinit <- matrix(runif(ncol(data)*3), ncol = ncol(data))
## Run traditional Tweedie NMF:
#res_init <- nmfgen(data, 3, "Tweedie", "traditional", pwr = 0, wmat = winit, hmat = hinit)
```

As before, the output can then be accessed as follows:

``` r
## Access weight matrix
#head(res_init$W)
## Access features matrix
#head(res_init$H)
```
