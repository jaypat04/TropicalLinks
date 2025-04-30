
function rand_nonzero(K)
    while true
        i = K(rand(-99:99))
        if i != zero(K)
            return i
        end
    end
end

function experiment(I::MPolyIdeal)
    R = base_ring(I) # Base ring
    indep_sets = Singular.independent_sets(Singular.std(Oscar.singular_generators(I)))
    indep_sets = [R.(indep_set) for indep_set in indep_sets if length(indep_set) == dim(I)] # Filter independent sets to only those of the correct dimension
        


    for indep_set in indep_sets
        set_score = score(I,indep_set) # Compute the score for the independent set

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

