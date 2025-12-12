
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CdaPrac1

<!-- badges: start -->

<!-- badges: end -->

CdaPrac1 is an educational R package designed to help students
understand and implement Coordinate Descent Algorithms (CDA),
particularly for the Lasso regression.

This package provides:

A pure R implementation of Lasso via Coordinate Descent

A fast Rcpp-based version (lasso_cda_cpp)

A user-friendly wrapper (lasso_cda) that performs automatic
standardization, coefficient recovery, and intercept calculation

Clean examples for learning and benchmarking

## Installation

You can install the development version of CdaPrac1 from GitHub:

``` r
# install.packages("pak")
pak::pak("LJWstat/CdaPrac1")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CdaPrac1)
# 
# set.seed(1)
# 
# # Simulated data
# 
# n <- 100
# p <- 10
# X <- matrix(rnorm(n*p), n, p)
# beta_true <- c(3, -2, rep(0, 8))
# y <- X %*% beta_true + rnorm(n)
# 
# # Fit Lasso using CDA (R + Rcpp implementation)
# 
# fit <- lasso_cda(X, y, lambda = 0.5)
# 
# fit$beta       # estimated coefficients
# fit$intercept  # estimated intercept
```

### Benchmark: R vs Rcpp

``` r
# library(microbenchmark)
# 
# fit_r <- function() lasso_cda_r(X, y, lambda = 0.5)
# fit_cpp <- function() lasso_cda_cpp(X, y, 0.5)
# 
# microbenchmark(
# R = fit_r(),
# Rcpp = fit_cpp(),
# times = 20
# )
```

### Plot Example

``` r
# grid <- seq(-1, 1, length.out = 200)
# 
# plot(grid, grid^2, type = "l",
# main = "Example Plot (replace with real Lasso visualization)",
# col = "steelblue", lwd = 2)
```
