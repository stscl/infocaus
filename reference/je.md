# Joint Entropy

Estimate the joint entropy of selected variables.

## Usage

``` r
je(data, indices, base = exp(1), type = c("cont", "disc"), k = 3)
```

## Arguments

- data:

  Observation data.

- indices:

  Integer vector of column indices to include in joint entropy
  calculation.

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
infoxtr::je(matrix(1:100,ncol=2),1:2)
#> [1] 5.784231
```
