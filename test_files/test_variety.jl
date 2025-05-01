using Oscar

# This is how to compute 0-dimensional Tropical Varieties
# There are faster ways, not official yet

R, (x,y) = polynomial_ring(QQ, ["x", "y"])

nu = tropical_semiring_map(QQ, 2) # p-adic valuation

f = x - 2
g = y - 2

I = ideal([f,g])

Vector{QQFieldElem}.(vertices(first(tropical_variety(I, nu)))) 
