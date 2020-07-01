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
E = compute_E(f, sym)
G = groebner(f, sym, order='grevlex')
print("This is E", E)
print(G)
print(G.is_zero_dimensional)
print("This is E in Q_0(f)", G.reduce(E.as_expr())[1])
E_red = G.reduce(E.as_expr())[1]
E_mon = E_red.terms()[0]
Bounds = [g.LM(order='grevlex').exponents for g in G]
Quotient = twovar_list(Bounds)
Q_Mon = [Monomial(exp, sym) for exp in Quotient]
print("dim of Q_0(f)", len(Q_Mon))
AA = []
for a in Q_Mon:
    Row = []
    for b in Q_Mon:
        c = a*b
        q = G.reduce(c.as_expr())[1]
        print("c", c, 'reduced', q)
        print(div(q, E_red)[0])
        print(phi(c, E_mon, G))
        Row.append(phi(c, E_mon, G))
    AA.append(Row)
AA = Matrix(AA)
print(AA)
#D = AA.eigenvals()
#print(D)
#print(partial_determinants(AA))
print(similarity_diagonalize(AA))

