# FoSIntro - convenience functions for function-on-scalar regression with pffr()

Overview
------------

This R package is the accompanying package to Bauer et al. (2017) - ''An introduction to semiparametric function-on-scalar regression''.

A vignette giving an overview of the functions of the package can be found <a href='http://bauer.userweb.mwn.de/files/FoSIntro.pdf' target='_blank'>here</a>.

All examples shown in the paper can be reproduced using the code in the `_examples` folder. The vignette instead does not reproduce each single plot of the paper, but instead gives an overview of the functions provided by the package. These are mainly convenient effect/prediction/residual plot functions for function-on-scalar models fitted with `refund::pffr`, but also comprise functions to perform (non)parametric bootstrapping to retrieve confidence intervals for smooth effects from `pffr` models. 


Setup
------------

Install from GitHub using:

``` r
devtools::install_github("bauer-alex/FoSIntro")
```
