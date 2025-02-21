# Function to make a vector primitive
function make_primitive(v::Vector{QQFieldElem})
    return v./ gcd([numerator.(v)]...) # Divides each element of the vector by the GCD of the elements
end

function tropical_link(I::MPolyIdeal , nu::TropicalSemiringMap)
    d = 1 # Assuming d = 1 for now
    R = base_ring(I) # Base ring
    n = ngens(R) # Number of variables

    W = Vector{QQFieldElem}[] # Initialise the result set

    for i in 1:n
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


