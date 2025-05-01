# 🌴 Tropical Links

This repository is part of my dissertation on computing tropical varieties. It includes an implementation of the **tropical link algorithm**, based on the methods described by [Hofmann and Ren (2018)](https://doi.org/10.1007/s00454-018-0023-z), using the [Oscar](https://oscar.computeralgebra.de) algebraic geometry system in Julia.

It also contains **two experiments** that explore heuristic methods for selecting key variable subsets in the tropical **link** and **point** algorithms proposed in the same paper.

---

## 📁 Files Overview

### `tropical_link_functions.jl`

This file contains the implementation of the tropical link algorithm as desccribed by [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z). It includes the following functions:

  - **`make_primitive`**  
    Converts a vector of rational numbers into its primitive form.
    
  - **`first_A`**  
    Computes the key variable subset required for the tropical link algorithm.

  - **`tropical_link`**  
    Computes the tropical link of an ideal, using the previously computed key variable subset.

---

### `tropical_link_experiment_functions.jl`

This file provides an experimental framework that explores a heuristic approach to selecting the key variable subset for the tropical link algorithm. It contains the following components:

  - `tropical_link_experiment`
    A modified version of the tropical link function that uses heuristics for the key variable subset.

  - **Heuristic scoring functions:**
    - `degree_count_score` - Scores variables based on degree counting.
    - `initial_score` - Uses the initial ideal.
  
  - **Condition functions:**
    - `is_A_done` - Checks if the key variable subset is complete.
    - `is_A_admissible` - Checks if the key variable subset is admissible for the link function.
  
  - `greedy_selection`
    The core function that uses the scoring and condition functions to greedily select the key variable subset. This selection is then used for computing the tropical link with the selected subset.

### `tropical_point_experiment_functions.jl`

This file focuses on experiment related to the tropical point algorithm in [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z) and explores a heuristic approach for choosing the maximal independent set module the ideal. It includes the following functions:

  - `rand_nonzero(K)`  
    Generates a random non-zero element in the specified field.

  - `score`
    A heuristic score function that scores an independent set using a degree counter.

  - `experiment`
    The main experiment function that:
    - Computes the maximal independent sets of an ideal
    - Scores the sets
    - Reduces the ideal to dimension 0 by mapping the variables in the independent set to a random element of the field.
    - Measures the time taken to compute the triangular decomposition of the 0-dimensional ideal.

  - `random_ideal`
    A function that generates a random ideal for testing the experiment.
    
  - `inflate_ideal(I::MPolyIdeal)`
    This function increases the exponents of certain randomly chosen variables in the ideal.

---

### `4. TropicalLinks.jl`

This file serves as the entry point to the project. It imports all required packages: Oscar, Combinatorics, and Random. It then loads the previously mentioned files including all the core and experimental functions. By running this script in Julia, all neccessary functions and depdencies are initialised, making it easy to access everything in one place.

---

## 🔬 References

[HR18] - Hofmann, T., & Ren, Y. (2018). Computing tropical points and tropical links. *Discrete & Computational Geometry, 60*(3), 627-645. https://doi.org/10.1007/s00454-018-0023-z
