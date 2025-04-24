R, (x1,x2,x3,x4,x5,x6,x7,x8,x9) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"]) # Define the polynomial ring
I = ideal([x1*x5 - x2*x6, x3*x4 - x1*x6, x7*x8 - x5*x9, x2*x8 - x3*x7])

A = first_A(I) # Get the first admissible set A
Acomplement = setdiff(1:ngens(R), A) # The variables to map to 1
scores = [degree_count_score(I,i) for i in Acomplement]
min_score_index = argmin(scores) # Get the index of the variable with the minimum score
deleteat!(Acomplement, min_score_index) # Remove the variable with the minimum score from Acomplement   

R = base_ring(I) # Base ring
keep_indices = setdiff(1:ngens(R), Acomplement) # The variables to keep
R0symbols = [copy(symbols(R))[i] for i in keep_indices]
R0, x0 = polynomial_ring(QQ, R0symbols) # Define the polynomial ring

phi = hom(R, R0, [i in Acomplement ? one(R0) : x0[findfirst(==(i), keep_indices)] for i in 1:ngens(R)]) # Define the homomorphism

J = phi(I) # Apply the homomorphism to the ideal
