# Returns a random non-zero element of the field K
function rand_nonzero(K)
    while true
        i = K(rand(-99:99))
        if i != zero(K)
            return i
        end
    end
end

# The score function computes the score of a given independent set of variables
function score(I::MPolyIdeal, vars::Vector{<:MPolyRingElem})
    score = 0

    for f in gens(I) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            for var in unique(vars) # Iterate over variables in the independent set               
                score += degree(term, var) # Increment the score by the degree of that term with respect to the variable
            end
        end
    end
    return score # Return the dictionary with variables and their scores
end

# Runs the experiment on the ideal I
function experiment(I::MPolyIdeal)
    R = base_ring(I) # Base ring
    indep_sets = Singular.independent_sets(Singular.std(Oscar.singular_generators(I))) # Get the independent sets from Singular
    indep_sets = [R.(indep_set) for indep_set in indep_sets if length(indep_set) == dim(I)] # Filter for maximal independent sets only
        


    for indep_set in indep_sets # Iterate over the independent sets
        set_score = score(I,indep_set) # Compute the score for the independent set

        # Constructing the homomorphism
        indep_set_indices = [findfirst(==(var), gens(R)) for var in indep_set] # The variables to map to 1
        keep_indices = setdiff(1:ngens(R), indep_set_indices) # The variables to keep
        reducedRsymbols = [copy(symbols(R))[i] for i in keep_indices]
        reducedR, reducedx = polynomial_ring(coefficient_ring(R), reducedRsymbols) # Define the polynomial ring    
        phi = hom(R, reducedR, [i in indep_set_indices ? rand_nonzero(coefficient_ring(R)) : reducedx[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism
    
        J = phi(I) # Apply the homomorphism to the ideal

        # Print the results
        println("The independent set: ", indep_set)
        println("Score = ", set_score)
        println("With time: ", @time triangular_decomposition(J))

    end
end


# This function generates a random ideal in a polynomial ring over QQ
# with a specified number of variables, generators, and maximum degree of the generators.
function random_ideal(rng=Random.default_rng(); nvars=12, n_gens=6, max_deg=3)
    R, _ = polynomial_ring(QQ, ["x$i" for i in 1:nvars]) # Create a polynomial ring with nvars variables

    # Generate random polynomials
    polys = []
    for _ in 1:n_gens
        monomials = gens(ideal(gens(R))^max_deg)
        push!(polys, sum(rand_nonzero(coefficient_ring(R)) * m for m in rand(monomials, 3)))
    end

    I = ideal([polys...]) # Form the ideal

    return(I)
end

# This function inflates the ideal by mapping a subset of variables to a higher power
function inflate_variable(I::MPolyIdeal)
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables
    k = ngens(I) # Number of generators

    random_indices = rand(collect(combinations(1:n, n-k))) # Randomly select n-k indices to inflate
    vars = gens(R)[random_indices] # The variables to map to a higher power

    phi = hom(R, R, [var in vars ? var^10 : var for var in gens(R)]) # Define the homomorphism
    J = phi(I) # Apply the homomorphism to the ideal
    
    return J 
end

