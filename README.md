# Tropical Links

This project is a part of my dissertation on computing tropical varieties.
It includes the implementation of the algorithm for computing tropical links, which is based on the methods discussed in [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z), in Oscar. 
As well as two experiments that explore possible ways of selecting key variable subsets in the tropical point and tropical link algorithm of [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z).

## Files Overview

### **1. tropical_link_functions.jl**

This file contains the implementation of the tropical link algorithm as desccribed by [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z). It includes the following functions:
  - **make_primitive(v::Vector{QQFieldElem})**: This function takes a vector of rational numbers and converts them to a primitive form.
    
  - **first_A(I::MPolyIdeal)**: This function computes the key variable subset required for the tropical link algorithm.

  - **tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap)**: This function computes the tropical link of an ideal, using the previously computed key variable subset.

### **2. tropical_link_experiment_functions.jl**

This file implements an experimental framework that explores a heuristic approach to selecting the key variable subset for the tropical link algorithm. It contains the following components:

  - **tropical_link_experiment(I::MPolyIdeal , nu::TropicalSemiringMap)**: A modified version of the tropical link function that uses heuristics for the key variable subset.

  - **Heuristic scoring functions**:
    - **degree_count_score(I::MPolyIdeal, index::Int64)**: A score algorithm that uses a degree counting method.
    - **initial_score(I::MPolyIdeal, nu::TropicalSemiringMap,  w::Vector, index::Int64)**: A similar score algorithm but uses the initial ideal instead.
  
  - **Condition functions**:
    - **is_A_done(I::MPolyIdeal, A::Vector{Int64})**: Checks if the key variable subset is complete.
    - **is_A_admissible(I::MPolyIdeal, A::Vector{Int64})**: Checks if the key variable subset is admissible for the link function.
  
  - **greedy_selection(I::MPolyIdeal, score::Function, is_admissible::Function, is_done::Function)**: The core function that uses the scoring and condition functions to greedily select the key variable subset. This selection is then used for computing the tropical link with the selected subset.





## References

[HR18] - Hofmann, T., & Ren, Y. (2018). Computing tropical points and tropical links. *Discrete & Computational Geometry, 60*(3), 627-645. https://doi.org/10.1007/s00454-018-0023-z
