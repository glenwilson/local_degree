from function_defs import *
from sympy import Expr, Poly, Symbol, groebner, itermonomials, Monomial, Matrix
from sympy.polys.orderings import monomial_key


n = 2
sym = [Symbol('x'+str(i)) for i in range(n)]
x = sym[0]
y = sym[1]
f = [Poly(-x**3*y+x*y**3), Poly(x**4+y**4-3*x**2*y**2)]
EKL = two_dim_EKL(f, sym)
print(EKL)
print("rank: ", len(EKL), "signature: ", signature(EKL))

for i in range(50):
    g = [rand_poly(5,4,1,sym), rand_poly(4,3,1,sym)]
    print(g)
    EKL = two_dim_EKL(g,sym)
    print(EKL)
    print("rank: ", len(EKL), "signature: ", signature(EKL))

