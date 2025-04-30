R9, (x1,x2,x3,x4,x5,x6,x7,x8,x9) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"]) # Define the polynomial ring
I = ideal([x1*x5 - x2*x6, x3*x4 - x1*x6, x7*x8 - x5*x9, x2*x8 - x3*x7])

R5, (x1,x2,x3,x4,x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"]) # Define the polynomial ring
I2 = ideal([2x1^3*x2 - 5x3^2*x4 + x5^3, x1^2*x3*x4 -3x2^4 + 4x5^2, 4x1*x2^2*x3 - 6x4^2*x5])
I3 = ideal([2x1^3*x2^2-5x2*x4^2+3x3^4-x5^5, 4x1^2*x3^3-2x2^5+x4^3*x5-3x1*x5^2, x3^4*x5 - x1^4*x4])
I4 = ideal([2x1^10*x2 - 5x3^6*x4 + x5^2, x1^2*x3*x4 -3x2^3 + 4x5^2, 4x1*x2^2*x3 - 6x4^5*x5])



R10, (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"]) # Define the polynomial ring

f1 = 2x1^4*x5^2-3x2^5*x7+x3^3*x6^4-5x4^6
f2 = -x2^3*x4^5+4x5^2*x6^3-x8^4+2x9^2*x10^3
f3 = x1^2*x7^5-3x4^3*x6^4+5x5^3-2x3^4*x9^2 + x10^4
I5 = ideal([f1,f2,f3])