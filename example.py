from sympy import *
x, y, z = symbols('x,y,z')
MyVars = [Symbol('x'+str(i)) for i in range(1,11)]
print(MyVars)
print(type(MyVars[0]))
e = MyVars[0]**2 - MyVars[1]**2
print(type(e))
print(e)
f = MyVars[0]*MyVars[1] + 2*MyVars[1]**2
e.as_poly(domain='QQ')
f.as_poly(domain='QQ')
print(type(e))
print(e)
#q, r = div(e, f, domain='QQ')
#print(q)
#print(r)
G=groebner([e,f], MyVars[0], MyVars[1], order='grevlex', domain='QQ')
print(G)
print(type(G))
print(G.is_zero_dimensional)
m = MyVars[0]**3 + MyVars[1]**2
m.as_poly(domain='QQ')
print("the type of m")
print(type(m))
print(G.reduce(m))
print(G._basis)
print(type(G._basis[0]))
print(G._basis[0])
print(type(G[0]))
print(G[0])

A = []
A.append(div(e,MyVars[0]))
A.append(div(A[0][1],MyVars[1]))
print(A)
