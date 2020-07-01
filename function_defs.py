from sympy import div, Matrix, itermonomials
from sympy.polys.orderings import monomial_key

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

def onevar_list(Bounds):
    n = min(Bounds)
    return range(n)

def twovar_list(Bounds):
    Bounds =sorted(Bounds, key=monomial_key(order='lex'))
    print(Bounds)
    out = []
    for i in range(len(Bounds)-1):
        b = Bounds[i]
        c = Bounds[i+1]
        out = out + [(x, y) for x in range(b[0], c[0]) for y in range(0, b[1])]
    return out

def region_list(Bounds):
    d = {}
    n = len(Bounds[0])
    for i in xrange(n):
        d[i] = range(max([b[i] for b in Bounds]))
    return output

def phi(test, E_monomial, G):
    """
    This ensures phi(E_monomial) = phi(E) = 1. 
    Returns phi(test)
    G is the groebner basis, 
    """
    test_red = G.reduce(test.as_expr())[1]
    for monomial in test_red.terms():
        if monomial[0] == E_monomial[0]:
            return monomial[1]*(E_monomial[1]**-1)
    return 0

def partial_determinants(M):
    """M is a square matrix, n x n, return a list of the determinants of
    the i x i principal minors for 1 <= i <= n.

    """
    n = M.rows
    factors = []
    if n != M.cols:
        raise TypeError("M is not square")
    for i in range(0,n):
        Sub = M.extract(range(0,i), range(0,i))
        factors.append(Sub.det())
    return factors

def similarity_diagonalize(M):
    """
    We assume M is symmetric and perform a diagonalization with simultaneous row and column operations.
    """
    n = M.rows
    diag = []
    if n != M.cols:
        raise TypeError("M is not square")
    #If M0,0 is non-zero, we can use it to clear the row and
    #column. Otherwise, we must find a non-zero entry in the first
    #row/col to put in this place.
    if M[(0,0)] == 0:
        for i in range(1,n):
            if M[(i,0)] != 0:
                #print('row ', i, ' has nonzero entry', M[(i,0)], M[(0,i)])
                M = M.elementary_row_op('n->n+km', k = 1, row=0, row2=i)
                M = M.elementary_col_op('n->n+km', k = 1, col=0, col2=i)
                break
    #Now M0,0 is non-zero so can clear the row and column
    #print('0,0-entry is ', M[(0,0)] )
    for i in range(1,n):
        if M[(i,0)] != 0:
            M = M.elementary_row_op('n->n+km', k = -M[(i,0)] *(M[(0,0)]**(-1)), row = i, row2 = 0)
            M = M.elementary_col_op('n->n+km', k = -M[(i,0)] *(M[(0,0)]**(-1)), col = i, col2 = 0)
        #print(M)
    #Now the first row and column are as desired. Repeat this procedure on the submatrix.
    diag.append(M[(0,0)])
    if n == 1:
        return diag
    else:
        return diag + similarity_diagonalize(M.extract(range(1,n),range(1,n)))
