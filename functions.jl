function tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap)
    A = variable_subset(I) # Get the variable subset
    R = base_ring(I) # Base ring
    keep_indices = setdiff(1:ngens(R), A) # Indices of the variables to keep
    R0symbols = [copy(symbols(R))[i] for i in keep_indices] # Define the symbols for the new polynomial ring
    R0, x0 = polynomial_ring(QQ, R0symbols) # Define the new polynomial ring
    phi = hom(R, R0, [i in A ? one(R0) : gens(R0)[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

    phi(I) # Apply the homomorphism to the ideal
    
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables

    W = Vector{QQFieldElem}[] # Initialise the result set

    for i in d:n
        p = base_ring(R)(uniformizer(nu)) # Uniformizer

        Risymbols = copy(symbols(R)) # Copy the list of symbols
        deleteat!(Risymbols, i) # Remove the i-th symbol

        Ri, xWithouti = polynomial_ring(QQ, Risymbols) # Define the polynomial ring without the i-th variable

        phiplus = hom(R, Ri, insert!(copy(xWithouti), i, Ri(p))) # Define the homomorphism for the positive map
        phiminus = hom(R, Ri, insert!(copy(xWithouti), i, Ri(p^(-1)))) # Define the homomorphism for the negative map

        Jplus = phiplus(I) # Apply the positive map
        Jminus = phiminus(I) # Apply the negative map

        Tplus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jplus, nu)))) # Compute the tropical variety for the positive map
        Tminus = Vector{QQFieldElem}.(vertices(first(tropical_variety(Jminus, nu)))) # Compute the tropical variety for the negative map

        Tplus = [make_primitive(insert!(vcat(zeros(QQ, d - 1), t), i, 1)) for t in Tplus] 
        Tminus = [make_primitive(insert!(vcat(zeros(QQ, d - 1), t), i, -1)) for t in Tminus]

        
        # Add to the result set
        append!(W, Tplus)
        append!(W, Tminus)
    end

    return unique(W) # Return the unique elements
end

function variable_subset(I::MPolyIdeal)
    G = collect(groebner_basis(I, complete_reduction = true))
    H = matrix(QQ, lineality_space(homogeneity_space(G)))
    R = echelon_form(H)
    Acomplement = [findfirst(!iszero, R[i, :]) for i in 1:nrows(R)]

    _,cols = size(R) # n
    all_indices = collect(1:cols)
    A = setdiff(all_indices, Acomplement)   
    return A
end

# Function to make a vector primitive
function make_primitive(v::Vector{QQFieldElem})
    return v./ gcd([numerator.(v)]...) # Divides each element of the vector by the GCD of the elements
end

function find_pivot_indices(R)
    pivot_indices = []
    for j in 1:size(R, 2)  # Iterate over columns
        for i in 1:size(R, 1)  # Iterate over rows
            if E[i, j] == 1 && all(R[k, j] == 0 for k in 1:size(R, 1) if k != i) # Check if the column has a leading 1
                push!(pivot_indices, j)
                break  # Move to the next column after finding a pivot
            end
        end
    end
    return pivot_indices
end

function degree_count_score(I::MPolyIdeal) # Old function
    R = base_ring(I)
    vars = symbols(R) # Extract the variables in the polynomial ring
    var_scores = Dict(var => 0 for var in vars) # Initialize a dictionary with variables as keys and scores as 0

    for f in gens(I) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            for (i,var) in enumerate(vars)
                var_scores[var] += degree(term, i) # Increment the score by the degree of that term with respect to the variable
            end
        end
    end
    var_scores = sort(collect(var_scores); by = x -> x[2])
    return var_scores # Return the dictionary with variables and their scores
end

function degree_count_score(I::MPolyIdeal, var::MPolyRingElem)
    score = 0
    for f in gens(I) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            score += degree(term, var) # Increment the score by the degree of that term with respect to the variable
        end
    end
    return  score # Return the dictionary with variables and their scores
end

function initial_score(I::MPolyIdeal, nu::TropicalSemiringMap,  w::Vector) # Old function
    initial_ideal = initial(I, nu, w)
    R = base_ring(initial_ideal)
    vars = symbols(R) # Extract the variables in the polynomial ring
    var_scores = Dict(var => 0 for var in vars) # Initialize a dictionary with variables as keys and scores as 0
    println(initial_ideal)

    for f in gens(initial_ideal) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            for (i,var) in enumerate(vars)
                var_scores[var] += degree(term, i) # Increment the score by the degree of that term with respect to the variable
            end
        end
    end 
    var_scores = sort(collect(var_scores); by = x -> x[2])
    return var_scores # Return the dictionary with variables and their scores
end

function initial_score(I::MPolyIdeal, nu::TropicalSemiringMap,  w::Vector, var::MPolyRingElem)
    R = base_ring(I) # Base ring
    index = findfirst(==(var), gens(R)) # Find the index of the variable in the polynomial ring
    initial_ideal = initial(I, nu, w) # Compute the initial ideal
    score = 0

    for f in gens(initial_ideal) # Iterate over the generators of the ideal
        for term in terms(f) # Iterate over each term in the generator
            score += degree(term, index) # Increment the score by the degree of that term with respect to the variable
        end
    end 
    return score # Return the dictionary with variables and their scores   
end

function is_A_done(I::MPolyIdeal, lambda::Vector{MPolyRingElem}) # Tropical link case
    R = base_ring(I) # Base ring
    n = nvars(R) # Number of variables
    d = dim(I) - 1 # Dimension of the ideal - 1
    if length(lambda) == n-d    
        return true
    else
        return false
    end    
end

score_closure = var -> degree_count_score(I,var)

# Greedy selection for the degree score
function test_greedy_selection(I::MPolyIdeal, score::Function, vars::Vector{QQMPolyRingElem}) #, is_admissible::Function, is_done::Function)
    min_score = inf
    best_var = nothing

    for var in vars
       println(score(I, var))
         if score(I, var) < min_score
              min_score = score(I, var)
              best_var = var
         end
    end

    return best_var # Return the variable with the minimum score
end

# Greedy selection for the initial score
function test_greedy_selection(I::MPolyIdeal, score::Function, nu::TropicalSemiringMap, w::Vector, vars::Vector{QQMPolyRingElem}) #, is_admissible::Function, is_done::Function)
    for var in vars
        println(score(I, nu, w, var))
     end
end

function greedy_selection(I::MPolyIdeal, score::Function, vars::Vector{QQMPolyRingElem}, is_admissible::Function, is_done::Function)
    lambda = Vector{QQMPolyRingElem}() # Initialize the lambda vector
    while !is_done(lambda, I) # While the lambda vector is not done
        min_score = Inf # Initialize the minimum score to infinity
        best_var = nothing

        for var in vars # Iterate over the variables
            current_score = score(I, var) # Get the score of the variable
            if current_score < min_score # If the score is less than the minimum score
                min_score = current_score # Update the minimum score
                best_var = var # Update the best variable
            end
        end

        push!(lambda, best_var) # Add it to the lambda vector
        if is_admissible(lambda) # If the variable is admissible
            deleteat!(vars, findfirst(==(best_var), vars)) # Remove it from the list of variables 
        else
            pop!(lambda) # Remove the last added variable from lambda
        end
    end
    return lambda # Return the lambda vector
end