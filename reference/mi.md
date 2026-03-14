# Mutual Information

Estimate the mutual information between target and interacting
variables.

## Usage

``` r
mi(
  data,
  target,
  interact,
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
infoxtr::mi(matrix(1:100,ncol=2),1,2)
#> [1] 2.979205
```
