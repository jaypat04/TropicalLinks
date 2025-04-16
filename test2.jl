#A = variable_subset(I3) # Get the variable subset
A = [1, 2] # Example variable subset
R = base_ring(I3) # Base ring
keep_indices = setdiff(1:ngens(R), A)
R0symbols = [copy(symbols(R))[i] for i in keep_indices]
R0, x0 = polynomial_ring(QQ, R0symbols) # Define the polynomial ring without the i-th variable

phi = hom(R, R0, [i in A ? one(R0) : gens(R0)[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

phi(I3)

for i in 1:ngens(R0)
    println("i: ", i)
end

for i in 3:ngens(R)
    println("i: ", i)
end
