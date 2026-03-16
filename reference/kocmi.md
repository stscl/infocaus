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
  seed = 123456789,
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

## References

Zhang, X., Chen, L., 2025. Quantifying interventional causality by
knockoff operation. Science Advances 11.

## Examples

``` r
popd_df = readr::read_csv(system.file("case/popd.csv",package = "spEDM"))
#> Rows: 2806 Columns: 7
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> dbl (7): lon, lat, popd, elev, tem, pre, slope
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
popd_df
#> # A tibble: 2,806 × 7
#>      lon   lat  popd  elev   tem   pre slope
#>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1  117.  30.5  780.     8  17.4 1528. 0.452
#>  2  117.  30.6  395.    48  17.2 1487. 0.842
#>  3  117.  30.8  261.    49  16.0 1456. 3.56 
#>  4  116.  30.1  258.    23  17.4 1555. 0.932
#>  5  116.  30.5  211.   101  16.3 1494. 3.34 
#>  6  117.  31.0  386.    10  16.6 1382. 1.65 
#>  7  117.  30.2  350.    23  17.5 1569. 0.346
#>  8  117.  30.7  470.    22  17.1 1493. 1.88 
#>  9  117.  30.6 1226.    11  17.4 1526. 0.208
#> 10  116.  30.9  137.   598  13.9 1458. 5.92 
#> # ℹ 2,796 more rows
# infoxtr::kocmi(popd_df, 3, 6, c(4,5,7))
```
