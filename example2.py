from function_defs import compute_E
from sympy import *

sym = [Symbol('x'+str(i)) for i in range(2)]
x1 = sym[0]
x2 = sym[1]

f1 = Poly(x1**2 - x2**2)
f2 = Poly(x1*x2)

print(f1)
print(f2)

A = compute_E([f1, f2], sym)

B = Matrix(A)
print(B)
print(type(B))
print(B.det())
