R, (x,y,u,v,w) = polynomial_ring(QQ, ["x", "y", "u", "v", "w"]) # Define the polynomial ring
I = ideal([x*u - y*v, x*w - u*v])

Acomplement = variable_subset(I)[1:end-1] # The variables to map to 1
R = base_ring(I) # Base ring
keep_indices = setdiff(1:ngens(R), Acomplement) # The variables to keep
R0symbols = [copy(symbols(R))[i] for i in keep_indices]
R0, x0 = polynomial_ring(QQ, R0symbols) # Define the polynomial ring

phi = hom(R, R0, [i in Acomplement ? one(R0) : x0[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

phi(I) # Apply the homomorphism to the ideal
