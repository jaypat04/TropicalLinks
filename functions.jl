# Function to make a vector primitive
function make_primitive(v::Vector{QQFieldElem})
    return v./ gcd([numerator.(v)]...) # Divides each element of the vector by the GCD of the elements
end

function tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap)
    d = 1 # Assuming d = 1 for now
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

function degree_count_score(I::MPolyIdeal)
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