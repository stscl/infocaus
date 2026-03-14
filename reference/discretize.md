# Discretization

Discretize a numeric vector into categorical classes using several
commonly used discretization methods. Missing values (`NA`/`NaN`) are
ignored and returned as class `0`.

## Usage

``` r
discretize(
  x,
  n = 5,
  method = "natural",
  large = 3000,
  prop = 0.15,
  seed = 123456789,
  thr = 0.4,
  iter = 100,
  bps = NULL,
  right_closed = TRUE
)
```

## Arguments

- x:

  A vector.

- n:

  (optional) Number of classes.

- method:

  (optional) Discretization method. One of `"sd"`, `"equal"`,
  `"geometric"`, `"quantile"`, `"manual"`, `"natural("jenks")"`, or
  `"headtail"("headtails")`.

- large:

  (optional) Threshold sample size for natural breaks sampling.

- prop:

  (optional) Sampling proportion used when `method = "natural"` and the
  input size exceeds `large`.

- seed:

  (optional) Random seed used for sampling in natural breaks.

- thr:

  (optional) Threshold used in the head/tail breaks algorithm.

- iter:

  (optional) Maximum number of iterations for head/tail breaks.

- bps:

  (optional) Numeric vector of manual breakpoints used when
  `method = "manual"`.

- right_closed:

  (optional) Logical. If `TRUE`, intervals are right-closed (e.g.,
  `(a, b]`). If `FALSE`, intervals are left-closed `[a, b)`.

## Value

A discretized integer vector.

## Note

If `x` is not numeric, it will be converted to integer categories via
[`as.factor()`](https://rdrr.io/r/base/factor.html).

## Examples

``` r
set.seed(42)
infoxtr::discretize(stats::rnorm(99,1,10), n = 5, method = 'natural')
#>  [1] 89 27 59 70 60 43 94 44 98 46 87 99  9 36 41 71 35  2  3 88 34  5 40 85 97
#> [26] 31 37  6 63 25 62 76 83 26 65  7 21 19  4 47 53 33 78 24 10 61 20 93 30 74
#> [51] 57 22 96 72 50 55 75 51  1 56 32 52 67 92 23 86 58 84 81 77 15 45 69 16 28
#> [76] 66 79 64 17 14 95 54 49 42 11 68 38 39 82 80 91 29 73 90 13 18 12  8 48
```
