# Sea-Turtle Migration Model

A stochastic, individual-count simulation of how migration from a cooler source beach affects turtle populations on a beach affected by warming.

File: `Sea-Turtle_Migration-Model(beehive).Rmd`

## Model question

How does the immigration rate from a cool source beach affect the number of usual (type `a`) and adapted (type `A`) turtles on a warming beach, and the total population size, after many generations?

## Turtle types

| Type | Description |
|-----------------------|-------------------------------------------------|
| `a` | Usual turtles. Arrive as immigrants from the cool beach. Their expected offspring number is reduced by `decay_rate`. |
| `A` | Adapted turtles. Arise only through mutation from `a`. Their growth depends on `sel_coeff - decay_rate` and is density-dependent. |

## Assumptions

- Warming is not modelled as temperature. `decay_rate` stands for its negative effect on reproduction.
- Migration is a single immigration rate. Routes and route-specific mortality are not modelled.
- All immigrants are type `a`. No adapted turtles come from the source beach.
- The source beach population (`N_source_a`) is constant.
- Generations do not overlap. Each generation fully replaces the previous one.
- Births, mutations and immigrants are Poisson-distributed; trimming to carrying capacity is binomial.

## Parameters

| Parameter | Meaning | Default |
|---------------------------|-----------------------|-----------------------|
| `N_init_a` | Starting number of type `a` turtles on the warming beach | 1000 |
| `N_init_A` | Starting number of type `A` turtles on the warming beach | 0 |
| `decay_rate` | Negative effect of warming on expected offspring | 0.1 |
| `sel_coeff` | Relative growth advantage of type `A` | 0.2 |
| `mut_rate` | Per-offspring rate at which `a` mutates to `A` | 0.001 |
| `mig_rate` | Fraction of the source population arriving each generation | 0 to 1 |
| `N_source_a` | Size of the source (cool beach) population | 1000 |
| `t_max` | Number of generations simulated | 1000 |
| `carrying_capacity` | Maximum population the warming beach supports | 1000 |
| `n_replicates` | Number of repeated runs per parameter set | 100 |

All parameter sets are stored in `parameters.csv`, one scenario per row that can be used to generate different scenario.

## Step by Step for Running the model !

### 1 Clean the workspace

rm(list=ls())

### 2 Set the directory 

setwd() #set the directory into the location where we put our data

### 3. Installing the packages

install.packages("beeswarm")

### 4. Load the packages 

library(beeswarm)

### 5 Load the data & parameters

## Order of events in one generation

1.  **Reproduction**
    - `a` offspring \~ Poisson(`N_a × (1 − decay_rate)`)
    - `A` offspring \~ Poisson(`N_A + N_A × (sel_coeff − decay_rate) × (1 − (N_a + N_A) / K)`)
2.  **Carrying capacity** – if total offspring \> `K`, offspring are reduced to `K` by binomial sampling.
3.  **Mutation** – `a → A` mutations \~ Poisson(`offsp_a × mut_rate`).
4.  **Immigration** – new `a` immigrants \~ Poisson(`N_source_a × mig_rate`).
5.  **Update and cap** – new `N_a` and `N_A` are calculated, and carrying capacity is applied again after immigration.

## Functions

| Function | What it does |
|------------------------------|------------------------------------------|
| `simulate_one_gen()` | Runs one generation and returns `c(N_a, N_A)` |
| `simulate_pop()` | Runs `t_max` generations and returns a matrix with columns `a` and `A`, rows = generations 0 to `t_max` |

## Outputs

1.  **Population dynamics plot** – total, `a` and `A` over generations for one run.
2.  **Migration effect plot** – final-generation `A`, `a` and total population across migration rates (0, 10⁻⁵ … 10⁰). Points show each replicate (beeswarm); lines show replicate means.
3.  **Frequency table** – total population and frequencies of `a` and `A` for every generation.

## Running the model from `parameters.csv`

``` r
params <- read.csv("parameters.csv")

# Single trajectory
p <- params[params$scenario == "single_run", ]
output <- simulate_pop(
  N_init_a = p$N_init_a, N_init_A = p$N_init_A,
  decay_rate = p$decay_rate, sel_coeff = p$sel_coeff,
  mut_rate = p$mut_rate, mig_rate = p$mig_rate,
  N_source_a = p$N_source_a, t_max = p$t_max,
  carrying_capacity = p$carrying_capacity
)

# Migration sweep
sweep <- params[grepl("^mig_", params$scenario), ]
migration_rates <- sweep$mig_rate
n_replicates   <- sweep$n_replicates[1]
```

## Requirements

- R (≥ 4.0)
- `beeswarm` package: `install.packages("beeswarm")`

## Limitations

- No explicit temperature, sex ratio, age structure or temperature-dependent sex determination.
- Type `A` can only arise by mutation, so adaptation is slow when `mut_rate` is low.
- Results are stochastic; use `set.seed()` for reproducible runs.
