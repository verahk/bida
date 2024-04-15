
# bida
This repository constitutes an R-package that implements in `R` the Bayesian IDA (BIDA) to categorical variables.
The BIDA method estimates intervention distributions and causal effects from observational data, under the assumption that the causal model is an unknown causal Bayesian network. 
More specifically, for each cause-effect pair in the considered system, a posterior distribution over the associated intervention distributions are computed by combining Bayesian estimation of intervention distributions through the so-called backdoor formula with Bayesian model averaging. 


## Installation

You can install the development version of bida from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("verahk/bida_cat")
```



## Examples
A toy_example that illustrate the functionality of the package is shown in `./inst/examples/toy_example.Rmd`. 
