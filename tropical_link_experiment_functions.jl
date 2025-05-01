# The implementation of the tropical_link function that uses the greedy selection
function tropical_link_experiment(I::MPolyIdeal , nu::TropicalSemiringMap)
    R = base_ring(I) # Base ring

    A = greedy_selection(I, degree_count_score, is_A_admissible, is_A_done) # Get the admissible set A
    Acomplement = setdiff(1:ngens(R), A) # The variables to map to 1
    scores = [degree_count_score(I,i) for i in Acomplement]
    min_score_index = argmin(scores) # Get the index of the variable with the minimum score
    deleteat!(Acomplement, min_score_index) # Remove the variable with the minimum score from Acomplement

    keep_indices = setdiff(1:ngens(R), Acomplement) # The variables to keep
    R0symbols = [copy(symbols(R))[i] for i in keep_indices]
    R0, x0 = polynomial_ring(QQ, R0symbols) # Define the new polynomial ring
    phi = hom(R, R0, [i in Acomplement ? one(R0) : x0[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism
    J = phi(I) # Apply the homomorphism to the ideal
    

    W = Vector{QQFieldElem}[] # Initialise the result set

    for i in 1:ngens(R0)
        p = base_ring(R0)(uniformizer(nu)) # Uniformizer

        Risymbols = copy(symbols(R0)) # Copy the list of symbols
        deleteat!(Risymbols, i) # Remove the i-th symbol

        Ri, xWithouti = polynomial_ring(QQ, Risymbols) # Define the polynomial ring without the i-th variable

        phiplus = hom(R0, Ri, insert!(copy(xWithouti), i, Ri(p))) # Define the homomorphism for the positive map
        phiminus = hom(R0, Ri, insert!(copy(xWithouti), i, Ri(p^(-1)))) # Define the homomorphism for the negative map

        Jplus = phiplus(I) # Apply the positive map
        Jminus = phiminus(I) # Apply the negative map

        Tplus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jplus, nu)))) # Compute the tropical variety for the positive map
        Tminus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jminus, nu)))) # Compute the tropical variety for the negative map

        Wplus = [make_primitive(insert!(vcat(zeros(QQ, d - 1), t), i, 1)) for t in Tplus] 
        Wminus = [make_primitive(insert!(vcat(zeros(QQ, d - 1), t), i, -1)) for t in Tminus]

        
        # Add to the result set
        append!(W, Wplus)
        append!(W, Wminus)
    end

    return unique(W) # Return the unique elements
end

# The scoring function that counts the degree for the variable at the given index
function degree_count_score(I::MPolyIdeal, index::Int64)
    score = 0
    for f in gens(I) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            score += degree(term, index) # Increment the score by the degree of that term with respect to the variable
        end
    end
    return  score # Return the dictionary with variables and their scores
end

# The second scoring function that uses the initial of the ideal
function initial_score(I::MPolyIdeal, nu::TropicalSemiringMap,  w::Vector, index::Int64)
    initial_ideal = initial(I, nu, w) # Compute the initial ideal
    score = 0    

    for f in gens(initial_ideal) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            score += degree(term, index) # Increment the score by the degree of that term with respect to the variable
        end
    end 
    return score # Return the dictionary with variables and their scores   
end

# Checks the size of the admissible set A
function is_A_done(I::MPolyIdeal, A::Vector{Int64}) 
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables
    d = dim(I) - 1 # Dimension of the ideal - 1
    if length(A) == n-d    
        return true
    else
        return false
    end    
end

# Checks if A is admissible
function is_A_admissible(I::MPolyIdeal, A::Vector{Int64})
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables
    G = collect(groebner_basis(I, complete_reduction = true)) 
    H = matrix(QQ, lineality_space(homogeneity_space(G))) # Compute the lineality space
    H_cone = cone(H) # Compute the cone of the lineality space

    L = matrix(QQ, [i == j ? QQ(1) : QQ(0) for j in A, i in 1:n]) # Create the matrix Lin(e_I : i in A)
    L_cone = cone(L) # Compute the cone of the matrix

    intersection = intersect(H_cone, L_cone) # Compute the intersection of the cones

    return dim(intersection) == 0 # Return true if the intersection is empty
end

# The greedy selection function
function greedy_selection(I::MPolyIdeal, score::Function, is_admissible::Function, is_done::Function)
    R = base_ring(I) # Base ring
    lambda = Vector{Int64}() # Initialize the lambda vector
    available_indices = collect(1:ngens(R)) # Create a list of available indices

    while !is_done(I, lambda) # While the lambda vector is not done
        min_score = Inf # Initialize the minimum score to infinity
        best_var_index = nothing

        for i in available_indices # Iterate over the variables
            current_score = score(I, i) # Get the score of the variable
            if current_score < min_score # If the score is less than the minimum score
                min_score = current_score # Update the minimum score
                best_var_index = i # Update the best variable
            end
        end

        push!(lambda, best_var_index) # Add it to the lambda vector
        if is_admissible(I, lambda) # If the variable is admissible
            deleteat!(available_indices, findfirst(==(best_var_index), available_indices)) # Remove it from the list of variables 
        else
            pop!(lambda) # Remove the last added variable from lambda
            deleteat!(available_indices, findfirst(==(best_var_index), available_indices)) # Remove it from the list of variables 
        end
    end
    return lambda # Return the lambda vector
end