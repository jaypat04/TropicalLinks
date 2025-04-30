for i in 1:3
    print("Iteration: ", i, "\n")
    I = random_ideal()
    J = inflate_variable(I)
    experiment(J)
    println("=====================================", "\n")
end