using Random

function random_positive_dim_ideal(rng=Random.default_rng(); nvars=6, n_gens=3, max_deg=6)
    # Create a polynomial ring
    R, vars = polynomial_ring(QQ, ["x$i" for i in 1:nvars])

    while true
        # Generate random polynomials
        polys = []
        for _ in 1:n_gens
            # Random degree <= max_deg
            deg = rand(rng, 1:max_deg)
            push!(polys, rand_polynomial(R, deg))
        end

        # Form the ideal
        I = ideal([polys...])

        # Check dimension
        if dim(I) > 0
            return I
        end
        # Otherwise retry
    end
end

# Helper: generate a random polynomial
function rand_polynomial(R, deg)
    mons = monomials(R, 0:deg)  # monomials up to degree deg
    coeffs = [rand(-5:5) for _ in mons]  # random integer coefficients
    return sum(c * m for (c, m) in zip(coeffs, mons) if c != 0)
end

function rand_polynomial(R, deg)
    n = ngens(R)
    mons = []
    for _ in 0:deg  # randomly generate 20 terms
        exponents = rand(0:deg, n)
        mon = prod(gens(R)[i]^exponents[i] for i in 1:n)
        push!(mons, mon)
    end
    coeffs = [rand(-5:5) for _ in mons]
    return sum(c * m for (c, m) in zip(coeffs, mons) if c != 0)
end

#R, (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"]) # Define the polynomial ring
#R, (x1,x2,x3,x4,x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"]) # Define the polynomial ring
