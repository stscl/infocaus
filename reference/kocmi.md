# KOCMI

Knockoff Conditional Mutual Information

## Usage

``` r
kocmi(
  data,
  target,
  agent,
  conds,
  knockoff,
  null_knockoff = NULL,
  type = c("cont", "disc"),
  nboots = 10000,
  k = 3,
  threads = 1,
  seed = 42,
  base = exp(1),
  method = "equal",
  contain_null = TRUE
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer vector of column indices for the target variables.

- agent:

  Integer vector of column indices for the source (agent) variables.

- conds:

  Integer vector of column indices for the conditioning variables.

- knockoff:

  Knockoff realizations constructed for the `agent` variable while
  keeping the `target` variable unchanged. Each column corresponds to
  one Monte Carlo knockoff sample generated using the remaining
  variables except the target.

- null_knockoff:

  (optional) Knockoff realizations generated under the null setting
  where all variables are jointly used to construct knockoffs. Each
  column represents one Monte Carlo sample. If `contain_null = FALSE`,
  this argument can be `NULL`.

- type:

  (optional) Estimation method: `"disc"` for discrete mutual information
  or `"cont"` for continuous mutual information (KSG estimator).

- nboots:

  (optional) Number of permutations used in the sign-flipping
  permutation test for evaluating the significance of the mean
  information difference.

- k:

  (optional) For `type = "cont"`, the number of nearest neighbors used
  by the continuous conditional mutual information estimator. For
  `type = "disc"`, the number of bins used for discretization.

- threads:

  (optional) Number of threads used.

- seed:

  (optional) Random seed used for permutation test.

- base:

  (optional) Logarithm base of the entropy. Defaults to `exp(1)` (nats).
  Use `2` for bits or `10` for dits.

- method:

  (optional) Discretization method. One of `"sd"`, `"equal"`,
  `"geometric"`, `"quantile"`, `"natural("jenks")"`, or
  `"headtail"("headtails")`.

- contain_null:

  (optional) Logical. If `TRUE`, the test statistic is computed using
  knockoffs generated under the null model (provided in
  `null_knockoff`). In this case the difference is defined as \\I(Y;
  X\_{null} \| Z) - I(Y; X\_{knockoff} \| Z)\\. If `FALSE`, the original
  conditional mutual information \\I(Y; X \| Z)\\ is used instead and
  compared against the knockoff estimates \\I(Y; X\_{knockoff} \| Z)\\.

## Value

A named numeric vector.

## Note

`kocmi` only support numeric data.

## References

Zhang, X., Chen, L., 2025. Quantifying interventional causality by
knockoff operation. Science Advances 11.

## Examples

``` r
set.seed(42)
kn1 = replicate(50, stats::rnorm(100))
kn2 = replicate(50, stats::rnorm(100))
mat = replicate(3, stats::rnorm(100))
infoxtr::kocmi(mat, 1, 2, 3, kn1, kn2)
#>     t_stat    p_value 
#> -0.2059844  0.1540000 
```
