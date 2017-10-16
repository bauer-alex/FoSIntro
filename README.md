# FoSIntro - convenience functions for function-on-scalar regression with pffr()

Overview
------------

This R package is the accompanying package to Bauer et al. (2018) - ''An introduction to semiparametric function-on-scalar regression''.

All examples from the package can be reproduced using the code in the `_examples` folder. The vignette instead does not reproduce each single plot of the paper, but instead gives an overview of the functions provided by the package. These are mainly convenient effect/prediction/residual plot functions for function-on-scalar models fitted with `refund::pffr`, but also comprise functions to perform (non)parametric bootstrapping to retrieve confidence intervals for smooth effects from `pffr` models. 


Setup
------------

Install from GitHub using:

``` r
# Note: Building the vignette might take several minutes
devtools::install_github("bauer-alex/FoSIntro", build_vignettes=TRUE)
```

Look at the vignette:
``` r
browseVignettes(package = "FoSIntro")
```
