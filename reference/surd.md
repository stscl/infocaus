# SURD

Synergistic-Unique-Redundant decomposition of causality

## Usage

``` r
# S4 method for class 'data.frame'
surd(
  data,
  target,
  agent,
  lag = 1,
  bin = 5,
  method = "equal",
  max.combs = 10,
  threads = 1,
  base = 2,
  normalize = TRUE
)

# S4 method for class 'sf'
surd(
  data,
  target,
  agent,
  lag = 1,
  bin = 5,
  method = "equal",
  max.combs = 10,
  threads = 1,
  base = 2,
  normalize = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
surd(
  data,
  target,
  agent,
  lag = 1,
  bin = 5,
  method = "equal",
  max.combs = 10,
  threads = 1,
  base = 2,
  normalize = TRUE
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer vector of column indices for the target variables.

- agent:

  Integer vector of column indices for the source (agent) variables.

- lag:

  (optional) Lag of the agent variables.

- bin:

  (optional) Number of discretization bins.

- method:

  (optional) Discretization method. One of `"sd"`, `"equal"`,
  `"geometric"`, `"quantile"`, `"natural("jenks")"`, or
  `"headtail"("headtails")`.

- max.combs:

  (optional) Maximum combination order.

- threads:

  (optional) Number of threads used.

- base:

  (optional) Logarithm base of the entropy. Defaults to `exp(1)` (nats).
  Use `2` for bits or `10` for dits.

- normalize:

  (optional) Logical; if `TRUE`, return normalized mutual information.

- nb:

  (optional) Neighbours list.

## Value

A list.

- vars:

  Character vector indicating the variable combination associated with
  each information component.

- types:

  Character vector indicating the information type of each component.

- values:

  Numeric vector giving the magnitude of each information component.

## Note

SURD only support numeric data.

## References

Martinez-Sanchez, A., Arranz, G. & Lozano-Duran, A. Decomposing
causality into its synergistic, unique, and redundant components. Nat
Commun 15, 9296 (2024).

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
infoxtr::surd(columbus, 1, 2:3)
#> $vars
#> [1] "V1"       "V2"       "V1_V2"    "V1_V2"    "InfoLeak"
#> 
#> $types
#> [1] "U"        "U"        "R"        "S"        "InfoLeak"
#> 
#> $values
#> [1] 0.06337975 0.06743535 0.46246679 0.40671810 0.58577187
#> 
```
