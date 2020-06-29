from function_defs import *
from sympy import Expr, Poly, Symbol, groebner

n = 4
sym = [Symbol('x'+str(i)) for i in range(n)]
f = [Poly(sym[i]**2) for i in range(n)]
E = compute_E(f, sym)
G = groebner(f, sym, order='grevlex')
print(E)
print(G)
print(G.is_zero_dimensional)
print(E, type(E))
print(G.reduce(E.as_expr()))
