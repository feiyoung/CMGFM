# CMGFM
High-Dimensional Covariate-Augmented Generalized Factor Model

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/CMGFM)](https://cran.r-project.org/package=CMGFM)
[![](https://cranlogs.r-pkg.org/badges/CMGFM?color=orange)](https://cran.r-project.org/package=CMGFM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/CMGFM?color=orange)](https://cran.r-project.org/package=CMGFM)
<!-- badges: end -->


Existing methods for multi-omics representation learning often lack interpretability or overlook critical omics-specific and additional information. To address these limitations and meet the practical demands, we introduce CMGFM, an interpretable multi-omics representation learning approach via covariate-augumented generalized factor model. CMGFM is designed to account for cross-modal heterogeneity, capture nonlinear dependencies among the data, incorporate additional information, and provide excellent interpretability while maintaining high computational efficiency.



Check out  our   [Package Website](https://feiyoung.github.io/CMGFM/index.html) for a more complete description of the methods and analyses. 

# Installation
"CMGFM" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.

```{Rmd}
## Method 1：
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/CMGFM")

## Method 2: install from CRAN
install.packages("CMGFM")

```



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Simulated data](https://feiyoung.github.io/CMGFM/articles/simu.html)

## Simulated codes
For the codes in simulation study, check the `simu_code` directory of the repo.


## News

CMGFM version 1.1 released! (2024-06-23) 


