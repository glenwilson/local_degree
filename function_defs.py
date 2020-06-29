from sympy import div, Matrix

def compute_E(f, sym):
    """

    This takes in a polynomial map f, which is a list of n polynomials
    in n variables. It returns a value of E in the polynomial ring in
    n variables. This can later be reduced in Q_0(f). 

    sym is a list of the symbols used as variables for f. (so has length n)
    
    """
    m = []
    n = len(sym)
    for i in range(n):
        A = []
        r = f[i]
        for j in range(n):
            d = div(r, sym[j])
            A.append(d[0])
            r = d[1]
        m.append(A)
    B = Matrix(m)
    return B.det()
    
