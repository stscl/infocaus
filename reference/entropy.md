# Entropy Estimation

Estimate the entropy of a vector using either category counts (for
discrete data) or a k-nearest neighbor estimator (for continuous data).

## Usage

``` r
entropy(vec, base = exp(1), type = c("cont", "disc"), k = 3)
```

## Arguments

- vec:

  A vector.

- base:

  (optional) Logarithm base of the entropy. Defaults to `exp(1)` (nats).
  Use `2` for bits or `10` for dits.

- type:

  (optional) Estimation method: `"disc"` for discrete entropy or
  `"cont"` for continuous entropy (KSG estimator).

- k:

  (optional) Number of nearest neighbors used by the continuous
  estimator. Ignored when `type = "disc"`.

## Value

A numerical value.

## Examples

``` r
set.seed(42)
infoxtr::entropy(stats::rnorm(100), type = "cont")
#> [1] 1.463547
infoxtr::entropy(sample(letters[1:5], 100, TRUE), base = 2, type = "disc")
#> [1] 2.298349
```
