using Random

function random_ideal(rng=Random.default_rng(); nvars=8, n_gens=4, max_deg=4)
    R, _ = polynomial_ring(QQ, ["x$i" for i in 1:nvars])

    # Generate random polynomials
    polys = []
    for _ in 1:n_gens
        push!(polys, rand_polynomial(R, max_deg))
    end

    # Form the ideal
    I = ideal([polys...])
    return(I)
end


function rand_polynomial(R, max_deg)
    monomials = gens(ideal(gens(R))^max_deg)
    return sum(rand_nonzero(coefficient_ring(R)) * m for m in rand(monomials, 3))
end

function inflate_variable(I::MPolyIdeal)
    R = base_ring(I) # Base ring
    n = ngens(R)
    k = ngens(I)

    random_indices = rand(collect(combinations(1:n, n-k)))
    vars = gens(R)[random_indices]

    println(vars)

    phi = hom(R, R, [var in vars ? var^15 : var for var in gens(R)])
    J = phi(I) # Apply the homomorphism to the ideal
    
    return J
end
