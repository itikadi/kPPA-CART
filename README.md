# kPPA-CART (kurtosis-based Projection Pursuit Analysis augmented with Classification and Regression Trees)

Here, we implement kurtosis-based projection pursuit analysis (kPPA) augmented with classification and regression trees (CART) as an approach to deal with -omics data exhibiting small effect sizes and low feature intensities. We specifically employ kPPA-CART to integrate and visualize data deriving from the most common -omics platforms: transcriptomics by RNA-sequencing (RNA-seq), epigenomics by DNA methylation chips or bisulfite sequencing, and proteomics by reverse-phase protein arrays (RPPA) or mass spectrometry (MS). This package is described in Bong et al. (2025); *Nucleic Acids Research*; https://doi.org/10.1093/nar/gkaf844.

## Installation

You can install the package directly from GitHub using `devtools`. First, ensure you have `devtools` installed by running:

```R
install.packages("devtools")
devtools::install_github("FabianBong/KPPACart")
```

## Getting Started 
The current **implementation** of kPPA-CART is in R (this repository).
The **output** from any of the implementations can be explored either in R or Python.

### Tutorial

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
