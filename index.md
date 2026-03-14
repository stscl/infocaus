# infoxtr

[![infoxtr website:
https://stscl.github.io/infoxtr/](reference/figures/infoxtr.png)](https://stscl.github.io/infoxtr/)

***Information**-Theoretic Measures for Revealing Variable
**Interactions***

*infoxtr* is an R package for analyzing variable interactions using
information-theoretic measures. Originally tailored for time series, its
methods extend seamlessly to spatial cross-sectional data. Powered by a
pure C++ engine with a lightweight R interface, the package also exposes
its headers for direct integration into other R packages.

> *Refer to the package documentation <https://stscl.github.io/infoxtr/>
> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=infoxtr) with:

``` r
install.packages("infoxtr", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/infoxtr) with:

``` r
install.packages("infoxtr",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/infoxtr)
  with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/infoxtr",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
