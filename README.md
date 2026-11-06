# Sea-Turtle Migration Model

A stochastic simulation of how immigration from a cooler source beach affects a sea-turtle population on a beach affected by warming.

**Question:** How does immigration from a source population affect the demographic recovery and the establishment of an adapted type in a local sea-turtle population exposed to warming-related stress?

## Project structure

| File | What it is for |
|--------------------|----------------------------------------------------|
| `Final Sea-Turtle_Migration-Model(beehive).Rmd` | Model code, plots, results table and a Shiny app. |
| `parameters.csv` | Named parameter sets (`baseline`, `high`, `medium`, `verylow`). |
| `README.md` | This file. |

## Requirements and installation

- R (≥ 4.0), RStudio recommended
- R packages `beeswarm`, `shiny` and `shinyWidgets`, installed once with:

``` r
install.packages(c("beeswarm", "shiny", "shinyWidgets"))
```

## How to run

1.  Put the `.Rmd` file and `parameters.csv` in the same folder.
2.  Open that folder as an RStudio Project, or set it as the working directory with `setwd("path/to/folder")`.
3.  Open the `.Rmd` file and click **Knit** to produce the plots and table.
4.  To use the Shiny app, run its chunk directly in RStudio (the green ▶ button). A normal Knit creates a static document, so the app will not run there.

The migration sweep runs 100 replicates × 1000 generations for each migration rate, so it can take a few minutes. Results are random; add `set.seed(42)` at the top for reproducible runs.

### Parameters and `parameters.csv`

`parameters.csv` holds named parameter sets in a column called `parameter_set`. Choose one in the `.Rmd`:

``` r
selected_set <- "high"   # baseline, high, medium or verylow
```

Note: the simulations currently use the values written directly in the `simulate_pop()` calls in the `.Rmd`. To change a scenario, edit those values.

## How the model works

**Turtle types**

- `a` – usual turtles, both local ones and immigrants. Warming reduces their offspring by `decay_rate`.
- `A` – adapted turtles. Arise only by mutation from `a`. Grow at `sel_coeff − decay_rate`, limited by space on the beach.

**One generation**

1.  Reproduction:
    - `a` offspring \~ Poisson(`N_a × (1 − decay_rate)`)
    - `A` offspring \~ Poisson(`N_A + N_A × (sel_coeff − decay_rate) × (1 − (N_a + N_A) / K)`)
2.  If offspring exceed `K`, they are reduced to `K` by binomial sampling.
3.  Mutation `a → A` \~ Poisson(`offsp_a × mut_rate`). Mutants move from `a` to `A`.
4.  Immigration of `a` \~ Poisson(`N_source_a × mig_rate`).
5.  Population is updated and capped at `K` again by binomial sampling.

**What to expect.** The simulations start with an empty beach (`N_init_a = 0`), so without immigration the population stays at zero. With immigrants, `a` levels off near `mig_rate × N_source_a / decay_rate` (e.g. `mig_rate = 0.1` fills the beach). This creates the key trade-off: more migration keeps the population up (*demographic rescue*), but immigrants also take up space and slow the spread of `A` (*swamping*).

**Parameters used in the main simulations**

| Parameter | Meaning | Value |
|-----------------------------|------------------------|-------------------|
| `N_init_a`, `N_init_A` | Starting numbers of `a` and `A` | 0, 0 |
| `decay_rate` | Negative effect of warming on offspring | 0.1 |
| `sel_coeff` | Growth advantage of `A` | 0.2 |
| `mut_rate` | Rate of `a → A` mutation per offspring | 0.001 |
| `mig_rate` | Fraction of the source population arriving each generation | 0 (single run); 10⁻⁵ to 1 (sweep) |
| `N_source_a` | Size of the source population (constant) | 1000 |
| `t_max` | Number of generations | 1000 |
| `carrying_capacity` (K) | Maximum population on the warming beach | 1000 |
| `n_replicates` | Runs per migration rate | 100 |

## Outputs

1.  **Results table** – numbers of `a` and `A` in the last six generations of the single run.
2.  **Population dynamics plot** – total, `a` and `A` over time for the single run.
3.  **Migration effect plot** – final `a`, `A` and total population across migration rates. Points are replicates (beeswarm); lines are means.
4.  **Shiny app** ("Evolutionary Rescue Simulation") – an interactive, simplified version of the model. Sliders set the starting population, decline rate, selective advantage, mutation probability and number of generations. It shows population trajectories, the share of populations that rebound, and the expected rescue probability (Orr & Unckless 2014).

## Assumptions and limitations

- Warming is not modelled as temperature; `decay_rate` stands for its effect on reproduction.
- All immigrants are type `a`; the source population is constant.
- Generations do not overlap. No age structure, sex ratio or temperature-dependent sex determination.
- `A` arises only by mutation, so adaptation is slow when `mut_rate` is low.
- One generation is roughly 20–30 turtle years, so 1000 generations is a very long time span.

## Authors and context

Created by Fadly Muhammad and Sara Costanzo for Theoretical Ecology and Evolution, University of Bern, Fall 2026.
