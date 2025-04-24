using Combinatorics



R = base_ring(I)
n = ngens(R)
d = dim(I) - 1
G = collect(groebner_basis(I, complete_reduction = true))
H = matrix(QQ, lineality_space(homogeneity_space(G)))
H_cone = cone(H)

candidate_As = collect(combinations(1:n, n - d))
count = 0

for A in candidate_As
    L = matrix(QQ, [i == j ? QQ(1) : QQ(0) for j in A, i in 1:n])
    L_cone = cone(L)
    intersection = intersect(H_cone, L_cone)

    if dim(intersection) == 0
        println("Admissible set A: ", A)
        global count
        count += 1
    else
        println("NOT admissible set A: ", A)
    end
end
println("Number of admissible sets A: ", count)

