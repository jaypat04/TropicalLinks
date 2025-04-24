R = base_ring(I)
n = ngens(R)
d = dim(I) - 1

candidate_As = collect(combinations(1:n, n - d))

println(candidate_As[findfirst(A -> is_A_admissible(I, A), candidate_As)]) # Finds the first admissible set in candidate_As

println(is_A_admissible(I, [1,2,5,6,7]))
count = 0

for A in candidate_As
    if is_A_admissible(I, A)
        println("Admissible set: ", A)
        global count += 1
    else
        println("Not admissible set: ", A)
    end
end
    
println("Total admissible sets: ", count)

#12457
function first_A(I::MPolyIdeal) # Function that just selects A one by one in order of the variables (No score algorithm)
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