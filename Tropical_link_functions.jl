# Function to make a vector primitive
function make_primitive(v::Vector{QQFieldElem})
    return v./ gcd([numerator.(v)]...) # Divides each element of the vector by the GCD of the elements
end

 # Function that just selects A one by one in order of the variables (No score algorithm)
function first_A(I::MPolyIdeal)
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables
    d = dim(I) - 1 # Dimension of the ideal - 1
    A = Vector{Int64}() # Initialize A as an empty vector

    for i in 1:ngens(R) # Iterate over the variables
        push!(A, i) # Add the variable to lambda
        if !is_A_admissible(I, A) # Check if the variable is admissible
            pop!(A) # Remove the last added variable from lambda
        end

        if length(A) == n - d # If the set A is complete
            break # Break the loop
        end
    end 
    return A # Return the set A
end


# The implementation of the tropical_link function
function tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap)
    R = base_ring(I) # Base ring

    A = first_A(I) # Get the admissible set A
    Acomplement = setdiff(1:ngens(R), A)[1:end-1] # The variables to map to 1    
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