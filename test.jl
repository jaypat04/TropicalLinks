include("TropicalLinks.jl")

R, (x,y) = polynomial_ring(QQ, ["x", "y"])

f = 1 + x + y + x^2*y + x*y^2 + x^2*y^2

I = ideal([f]) # Input

nu = tropical_semiring_map(QQ, 3) # p-adic valuation

tropical_link(I, nu)