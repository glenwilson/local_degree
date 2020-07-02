from function_defs import *
from sympy import Expr, Poly, Symbol, groebner, itermonomials, Monomial, Matrix
from sympy.polys.orderings import monomial_key
from random import randint


n = 2
sym = [Symbol('x'+str(i)) for i in range(n)]
x = sym[0]
y = sym[1]
f = [Poly(-x**3*y+x*y**3), Poly(x**4+y**4-3*x**2*y**2)]
EKL = two_dim_EKL(f, sym)
print(EKL)
print("rank: ", len(EKL), "signature: ", signature(EKL))

# h = [f[0]*f[1], -f[0]**2+f[1]**2]
# EKL = two_dim_EKL(h, sym)
# print(EKL)
# print("rank: ", len(EKL), "signature: ", signature(EKL))


for i in range(1000):
    deg = 4
    P = Poly('0', sym)
    for i in range(deg+1):
        P += randint(-10,10) * Poly(x**(deg - i) * y**i)
#    deg = 3
#    for i in range(deg +1):
#        P += randint(-10,10) * Poly(x**(deg - i) * y**i)
    Q = Poly('0', sym)
    deg = 4
    for i in range(deg +1):
        Q += randint(-10,10) * Poly(x**(deg - i) * y**i)
    g = [Q, P]

    try:
        EKL = two_dim_EKL(g,sym)
        print(P,Q)
        print(EKL)
        print("rank: ", len(EKL), "signature: ", signature(EKL))
    except TypeError:
        continue
