from sympy import *
x, y, z = symbols('x,y,z')

e = x**2 - y**2
f = x*y
e.as_poly(domain='QQ')
f.as_poly(domain='QQ')
print(e)
print(f)
q, r = div(e, f, domain='QQ')
print(q)
print(r)
G=groebner([e,f], x, y, order='grevlex', domain='QQ')
print(G)
