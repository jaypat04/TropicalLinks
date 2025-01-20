include("TropicalLinks.jl")

R, (x,y) = polynomial_ring(QQ, ["x", "y"])

f = 1 + x + y + x^2*y + x*y^2 + x^2*y^2

I = ideal([f]) # Input

nu = tropical_semiring_map(QQ, 3) # p-adic valuation

tropical_link(I, nu)

i = 2
p = base_ring(R)(uniformizer(nu))

Risymbols = copy(symbols(R)) # Copy the list of symbols
deleteat!(Risymbols, i) # Remove the i-th symbol

Ri, xWithouti = polynomial_ring(QQ, Risymbols)

phiplus = hom(R, Ri, insert!(copy(xWithouti), i, Ri(p)))

phiminus = hom(R, Ri, insert!(copy(xWithouti), i, Ri(p^(-1))))

Jplus = phiplus(I)
Jminus = phiminus(I)

Tplus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jplus, nu)))) 
Tminus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jminus, nu))) )

