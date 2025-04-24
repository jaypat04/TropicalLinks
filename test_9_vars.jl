R, (x1,x2,x3,x4,x5,x6,x7,x8,x9) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"]) # Define the polynomial ring
I = ideal([x1*x5 - x2*x6, x3*x4 - x1*x6, x7*x8 - x5*x9, x2*x8 - x3*x7])

Acomplement = variable_subset(I)[1:end-1] # The variables to map to 1
R = base_ring(I) # Base ring
keep_indices = setdiff(1:ngens(R), Acomplement) # The variables to keep
R0symbols = [copy(symbols(R))[i] for i in keep_indices]
R0, x0 = polynomial_ring(QQ, R0symbols) # Define the polynomial ring

phi = hom(R, R0, [i in Acomplement ? one(R0) : x0[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

J = phi(I) # Apply the homomorphism to the ideal
