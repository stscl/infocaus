# Conditional Mutual Information

Estimate the conditional mutual information between target and
interacting variables given conditioning variables.

## Usage

``` r
cmi(
  data,
  target,
  interact,
  conds,
  base = exp(1),
  type = c("cont", "disc"),
  k = 3,
  normalize = FALSE
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer vector of column indices for the target variables.

- interact:

  Integer vector of column indices for the interacting variables.

- conds:

  Integer vector of column indices for the conditioning variables.

- base:

  (optional) Logarithm base of the entropy. Defaults to `exp(1)` (nats).
  Use `2` for bits or `10` for dits.

- type:

  (optional) Estimation method: `"disc"` for discrete entropy or
  `"cont"` for continuous entropy (KSG estimator).

- k:

  (optional) Number of nearest neighbors used by the continuous
  estimator. Ignored when `type = "disc"`.

- normalize:

  (optional) Logical; if `TRUE`, return normalized mutual information.

## Value

A numerical value.

## Examples

``` r
set.seed(42)
infoxtr::cmi(matrix(stats::rnorm(99,1,10),ncol=3),1,2,3)
#> [1] 0.02440107
```
