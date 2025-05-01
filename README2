# 🌴 Tropical Links

This repository is part of my dissertation on computing tropical varieties. It includes an implementation of the **tropical link algorithm**, based on the methods described by [Hofmann and Ren (2018)](https://doi.org/10.1007/s00454-018-0023-z), using the [Oscar](https://oscar.computeralgebra.de) algebraic geometry system in Julia.

It also contains **two experiments** that explore heuristic methods for selecting key variable subsets in the tropical **link** and **point** algorithms proposed in the same paper.

---

## 📁 Files Overview

### `tropical_link_functions.jl`

Implements the tropical link algorithm described in [Hofmann and Ren (2018)](https://doi.org/10.1007/s00454-018-0023-z). This file includes:

- **`make_primitive(v::Vector{QQFieldElem})`**  
  Converts a vector of rational numbers into its primitive form.

- **`first_A(I::MPolyIdeal)`**  
  Computes the key variable subset required by the tropical link algorithm.

- **`tropical_link(I::MPolyIdeal, ν::TropicalSemiringMap)`**  
  Computes the tropical link of an ideal using the previously computed key variable subset.

---

### `tropical_link_experiment_functions.jl`

Provides an experimental framework to explore **heuristic methods** for selecting the key variable subset used in the tropical link algorithm.

- **`tropical_link_experiment(I::MPolyIdeal, ν::TropicalSemiringMap)`**  
  A modified version of the tropical link algorithm that incorporates heuristics for choosing the variable subset.

- **Heuristic scoring functions:**
  - `degree_count_score(I::MPolyIdeal, index::Int64)` – Scores variables based on degree counting.
  - `initial_score(I::MPolyIdeal, ν::TropicalSemiringMap, w::Vector, index::Int64)` – Scores variables based on the initial ideal.

- **Condition functions:**
  - `is_A_done(I::MPolyIdeal, A::Vector{Int64})` – Checks if a key variable subset is complete.
  - `is_A_admissible(I::MPolyIdeal, A::Vector{Int64})` – Checks if the subset is admissible.

- **`greedy_selection(I::MPolyIdeal, score::Function, is_admissible::Function, is_done::Function)`**  
  A greedy algorithm that selects a key variable subset using scoring and condition functions.

---

### `tropical_point_experiment_functions.jl`

Focuses on a heuristic method for selecting a **maximal independent set modulo an ideal**, as part of the tropical point algorithm in [Hofmann and Ren (2018)](https://doi.org/10.1007/s00454-018-0023-z).

- **`rand_nonzero(K)`**  
  Generates a random non-zero element from a given field.

- **`score(I::MPolyIdeal, vars::Vector{<:MPolyRingElem})`**  
  Scores independent sets using a degree-based heuristic.

- **`experiment(I::MPolyIdeal)`**  
  The main experiment function:
  - Computes maximal independent sets modulo the ideal.
  - Scores each set.
  - Reduces the ideal to dimension 0 via random substitutions.
  - Times how long it takes to compute the triangular decomposition of the resulting ideal.

- **`random_ideal(rng=Random.default_rng(); nvars=12, n_gens=6, max_deg=3)`**  
  Generates a random ideal for testing.

- **`inflate_ideal(I::MPolyIdeal)`**  
  Increases the exponents of randomly chosen variables to modify the ideal structure.

---

### `TropicalLinks.jl`

Acts as the **entry point** for the project. This script:

- Imports required dependencies: `Oscar`, `Combinatorics`, and `Random`.
- Loads all source files:
  - `functions.jl`
  - `tropical_link_functions.jl`
  - `tropical_link_experiment_functions.jl`
  - `tropical_point_experiment_functions.jl`

By running this script, all algorithms and experiments are ready to use in a single Julia session.

---

## 🔬 Reference

> **[HR18]** T. Hofmann and Y. Ren. “Computing tropical points and tropical links.”  
> *Discrete & Computational Geometry*, 60(3), pp. 627–645, 2018.  
> [https://doi.org/10.1007/s00454-018-0023-z](https://doi.org/10.1007/s00454-018-0023-z)

