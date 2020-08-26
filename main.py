from function_defs import *
from sympy import Expr, Poly, Symbol, groebner, itermonomials, Monomial, Matrix
from sympy.polys.orderings import monomial_key
from random import randint


# n = 3
# sym = [Symbol('x'+str(i)) for i in range(n)]
# x = sym[0]
# y = sym[1]
# z = sym[2]
# f = [Poly(x**2+y**2+z**2), Poly(2*x*y*z), Poly(y**3+z**2)]
# EKL = three_dim_EKL(f, sym)
# print(EKL)
# print("rank: ", len(EKL), "signature: ", signature(EKL))

# ####################

# for i in range(1):
#     P = rand_poly(3,2,2,sym) 
#     Q = rand_poly(3,2,2,sym) 
#     R = rand_poly(3,2,2,sym) 
#     try:
#         EKL = three_dim_EKL([P,Q,R], sym)
#         print([P,Q,R])
#         print("rank of EKL", len(EKL))
#         print("signature of EKL", signature(EKL))
#         print("EKL", EKL)
#     except TypeError:
#         print([P,Q,R])
#         continue

####################

n = 2
sym = [Symbol('x'+str(i)) for i in range(n)]
x = sym[0]
y = sym[1]

f = [Poly(x**3-x**2+y**2), Poly((x-y)*(x+2*y))]
EKL = two_dim_EKL(f, sym)
print(EKL)
print("rank: ", len(EKL), "signature: ", signature(EKL))

for i in range(1000):
    s = rand_poly(1,1,5,sym)
    t = rand_poly(1,1,5,sym)
    u = rand_poly(1,1,5,sym)
    P = rand_poly(4,3,5,sym) + s*t
    Q = rand_poly(4,3,5,sym) + s*u
    try:
        EKL = two_dim_EKL([P,Q], sym)
        #print([P,Q])
        print("rank of EKL", len(EKL), "signature", signature(EKL))
        #print("signature of EKL", signature(EKL))
        #print("EKL", EKL)
    except TypeError:
        print("failure", [P,Q])
        continue
