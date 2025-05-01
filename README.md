# Tropical Links

This project is a part of my dissertation on computing tropical varieties.
It includes the implementation of the algorithm for computing tropical links, which is based on the methods discussed in [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z), in Oscar. 
As well as two experiments that explore possible ways of selecting key variable subsets in the tropical point and tropical link algorithm of [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z).

## Files Overview

### 1.tropical_link_functions.jl

This file contains the implementation of the tropical link algorithm as desccribed by [Hofmann and Ren, 2018](https://doi.org/10.1007/s00454-018-0023-z). It includes the following functions:
  - make_primitive(v::Vector{QQFieldElem}): This function takes a vector of rational numbers and converts them to a primitive form.
    
  - function first_A(I::MPolyIdeal): This function computes the key variable subset required for the tropical link algorithm.

  - tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap): This function computes the tropical link of an ideal, using the previously computed key variable subset.




## References

[HR18] - Hofmann, T., & Ren, Y. (2018). Computing tropical points and tropical links. *Discrete & Computational Geometry, 60*(3), 627-645. https://doi.org/10.1007/s00454-018-0023-z
