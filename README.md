# Iterative KPPA

This package adds iterative kPPA as described in "...".

## Installation

You can install the package directly from GitHub using the `devtools` package. First, ensure you have `devtools` installed by running:

```R
install.packages("devtools")
devtools::install_github("FabianBong/kPPA")
```

## Getting Started 

To see the functionality of the function it is easiest to load the sample data.

```R
data <- kPPA_sample_data()
res <- master_kPPA(data$X)
```

It is most helpful to plot the resulting scores in a 2/3-D plot.

```R
plot(res$T[,c(1,2)], col=data$Age)
```
