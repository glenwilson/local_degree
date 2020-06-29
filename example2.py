from function_defs import compute_E
from sympy import Poly, Symbol, groebner

sym = [Symbol('x'+str(i)) for i in range(3)]
x1 = sym[0]
x2 = sym[1]
x3 = sym[2]

f1 = Poly(x1**2 - x2**2)
f2 = Poly(x1*x2)
f3 = Poly(x3**3 + 2*x2**2)
print(f1, f2, f3)

E = compute_E([f1, f2, f3], sym)
m = Poly(x3**3, sym)
print("this is E", E)
G = groebner([f1, f2, f3], sym, order='grevlex')
print(G)
print(type(E.as_expr(*sym)), E)
print(G.reduce(E.as_expr()))
print(G.reduce(m.as_expr()))
