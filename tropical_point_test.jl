R, (x1,x2,x3,x4,x5,x6,x7,x8,x9) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"]) # Define the polynomial ring
I = ideal([x1*x5 - x2*x6, x3*x4 - x1*x6, x7*x8 - x5*x9, x2*x8 - x3*x7])

indep_sets = Singular.independent_sets(Singular.std(Oscar.singular_generators(I)))
indep_sets = [R.(indep_set) for indep_set in indep_sets]


function rand_nonzero_qq()
    r = rand(QQ, -10:10)
    while r == 0
        r = rand(QQ, -10:10)
    end
    return r
end


R = base_ring(I) # Base ring
indep_set_indices = [findfirst(==(var), gens(R)) for var in indep_sets[1]] # The variables to map to 1
keep_indices = setdiff(1:ngens(R), indep_set_indices) # The variables to keep


reducedRsymbols = [copy(symbols(R))[i] for i in keep_indices]
reducedR, reducedx = polynomial_ring(QQ, reducedRsymbols) # Define the polynomial ring

phi = hom(R, reducedR, [i in indep_set_indices ? rand_nonzero_qq() : reducedx[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

J = phi(I) # Apply the homomorphism to the ideal


