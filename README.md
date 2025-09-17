# KPPACart

This package adds iterative kPPA with classification and regression trees as described in Bong et al. (2025); Nucleic Acids Research; https://doi.org/10.1093/nar/gkaf844.

## Installation

You can install the package directly from GitHub using the `devtools` package. First, ensure you have `devtools` installed by running:

```R
install.packages("devtools")
devtools::install_github("FabianBong/KPPACart")
```

## Getting Started 

To see the functionality of the function it is easiest to load the sample data.

```R
data(KPPACart.Data)
res <- KPPACart(KPPACart.Data$X)
```

It is most helpful to plot the resulting scores in a 2/3-D plot.

```R
plot(res$T[,c(1,2)], col=KPPACart.Data$Age)
```
