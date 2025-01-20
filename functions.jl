function foo()
    println("Hello, World!")
end



#Assuming d = 1
function tropical_link(I, nu)
    R = base_ring(I)
    n = ngens(R) # Number of variables

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
        
    end

end
