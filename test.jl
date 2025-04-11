include("TropicalLinks.jl")

R, (x,y) = polynomial_ring(QQ, ["x", "y"])
f = 1 + x + y + x^2*y + x*y^2 + x^2*y^2

I = ideal([f]) # Input

nu = tropical_semiring_map(QQ, 3) # p-adic valuation

# tropical_link(I, nu)

W = []

d = 1
i =  2   
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

Tplus = [insert!(vcat(zeros(d - 1), t), i, 1) for t in Tplus] 
Tminus = [insert!(vcat(zeros(d - 1), t), i, -1) for t in Tminus]

append!(W, Tplus)
append!(W, Tminus)

R2, (x1,x2,x3) = polynomial_ring(QQ, [:x1, :x2, :x3])

# G1 = gens(I) # Generators of the ideal 

G = collect(groebner_basis(I, complete_reduction = true))

Gtest = [x1^2 + x2]
H = matrix(QQ, lineality_space(homogeneity_space(Gtest)))
R = echelon_form(H)

# E, pivots = echelon_form_with_transformation(H)  # Get echelon form and pivot indices

Acomplement = [findfirst(!iszero, R[i, :]) for i in 1:nrows(R)]

_,cols = size(R) # n
all_indices = collect(1:cols)
A = setdiff(all_indices, Acomplement)


R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"])

f = x^2*y + x^100*y
g = x + y^2 + x^3 + z

I = ideal([f, g]) # Input

vars = symbols(R) # Extract the variables in the polynomial ring
var_scores = Dict(var => 0 for var in vars) # Initialize a dictionary with variables as keys and scores as 0

for f in gens(I) # Iterate over the generators of the ideal
    for term in terms(f) # Iterate over each term in the generator
        for (i,var) in enumerate(vars)
            term_degree = degree(term, i) # Get the degree of the term with respect to the variable
            var_scores[var] += term_degree # Increment the score by the degree of that term
        end
    end
end

var_scores = sort(collect(var_scores); by = var -> var[2])
sorted_vars_list = [kv[1] for kv in var_scores]

println(var_scores) # Print the scores for each variable
println(sorted_vars_list) # Print the sorted list of variables based on their scores