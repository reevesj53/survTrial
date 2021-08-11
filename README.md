
<!-- README.md is generated from README.Rmd. Please edit that file -->

# survTrial

<!-- badges: start -->
<!-- badges: end -->

The goal of survTrial is to …

## Installation

You can install the released version of survTrial from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("reevesj53/survTrial")
```

## Example

This is a basic example which shows you how to solve a common problem:

Next, run simulation:

``` r
enrol <- c(seq(2,10,length.out=5),rep(10,times=3))
schedule <- seq(0,100,4)
rxrate <- c(12,10)
nevent <- 40
```
