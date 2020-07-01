from function_defs import *
from sympy import Expr, Poly, Symbol, groebner, itermonomials, Monomial, Matrix
from sympy.polys.orderings import monomial_key

n = 2
sym = [Symbol('x'+str(i)) for i in range(n)]
x = sym[0]
y = sym[1]
#f = [Poly(sym[i]**2) for i in range(n)]
#f[0] = f[0] + Poly(-sym[1]*sym[0])
f = [Poly(-(x**3)*y + x * y**3), Poly(x**4+y**4-3*x**2*y**2)]
print(two_dim_EKL(f, sym))
