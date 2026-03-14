# Discretization

Discretizes a numeric vector into categorical classes using several
commonly used discretization methods. Missing values (`NA`/`NaN`) are
ignored when computing breakpoints and returned as class `0`.

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

If `x` is a character or factor vector, it will be converted to integer
categories via [`as.factor()`](https://rdrr.io/r/base/factor.html).

## Examples

``` r
infocaus::discretize(1:100, n = 5, method = 'natural')
#> Error in loadNamespace(x): there is no package called ‘infocaus’
```
