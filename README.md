# KPPACart

This package adds iterative kPPA with classification and regression trees as described in Bong et al. (2025); *Nucleic Acids Research*; https://doi.org/10.1093/nar/gkaf844.

## Installation

You can install the package directly from GitHub using the `devtools` package. First, ensure you have `devtools` installed by running:

```R
install.packages("devtools")
devtools::install_github("FabianBong/KPPACart")
```

## Getting Started 
The current **implementation** of kPPA-CART is in R (this repository).
The **output** from any of the implementations can be explored either in R or Python.

### Tutorials

To use **R** for the whole analysis, there is the main tutorial:
> The tutorials in R include a more detailed explanation of the workflow and source code.

  - [Introduction and setup](https://htmlpreview.github.io/?https://github.com/itikadi/kPPA-CART/blob/main/vignettes/kPPA-CART.nb.html)

To see the functionality of the function, it is easiest to load the sample data.

```R
data(KPPACart.Data)
res <- KPPACart(KPPACart.Data$X)
```

It is most helpful to plot the resulting scores in a 2/3-D plot.

```R
plot(res$T[,c(1,2)], col=KPPACart.Data$Age)
```
