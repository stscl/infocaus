# Transfer Entropy

Estimate the transfer entropy from agent variables to target variables.

## Usage

``` r
te(
  data,
  target,
  agent,
  lag_p = 3,
  lag_q = 3,
  base = exp(1),
  type = c("cont", "disc"),
  k = 3,
  normalize = FALSE,
  lag_single = FALSE
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer vector of column indices for the target variables.

- agent:

  Integer vector of column indices for the source (agent) variables.

- lag_p:

  (optional) Lag of the target variables.

- lag_q:

  (optional) Lag of the agent variables.

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

- lag_single:

  (optional) Logical; if `FALSE`, use full lag embedding.

## Value

A numerical value.

## Examples

``` r
set.seed(42)
infoxtr::te(matrix(stats::rnorm(100,1,10),ncol=2),1,2)
#> [1] 0.01700429
```
