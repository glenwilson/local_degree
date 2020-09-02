from sympy import div, Matrix, itermonomials, groebner, Poly, Monomial, Expr
from sympy.polys.orderings import monomial_key
from random import randint
from sympy.polys.monomials import itermonomials

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
    out = []
    for i in range(len(Bounds)-1):
        b = Bounds[i]
        c = Bounds[i+1]
        out = out + [(x, y) for x in range(b[0], c[0]) for y in range(0, b[1])]
    return out

def threevar_list(Bounds):
    Bounds = sorted(Bounds, key=monomial_key(order='lex'))
    Max = max([a[0] for a in Bounds])
    #print(Max, "max")
    out = []
    for i in range(Max):
        sub_bounds = [a[1:3] for a in Bounds if a[0] <= i]
        #print(sub_bounds)
        #print(twovar_list(Bounds))
        out = out + [(i,) + x for x in twovar_list(sub_bounds)]
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


def groebner_m_test(n, G, sym):
    L = itermonomials(sym, n, 0)
    for l in L:
        test = Poly(0, sym)
        if (test + l).total_degree() == n and G.reduce(l)[1] != 0:
            print(l)
            print(G.reduce(l)[1])

def two_dim_EKL(f, sym):
    """Input is a polynomial function f, which in this case is a list of
    two polynomials in two variables. Sym is a list of the symbols
    used as variables in the polynomial.

    The output is the diagonal entries in the EKL class
    """
    E = compute_E(f, sym)
    G = groebner(f, sym, order='grevlex')
    if not G.is_zero_dimensional:
        raise TypeError("The polynomial chosen does not yield a 0-dimensional variety")
    #Bounds tells us the outer limits of the basis for the quotient 
    Bounds = [g.LM(order='grevlex').exponents for g in G]
    Quotient = twovar_list(Bounds)
    N = len(Quotient)
    #print(N, G)
    #print(Quotient)
    #N is the dimension of the quotient ring as a k-vector space. This
    #enables us to sort out the localization at 0 by further modding
    #out.
    f_loc = f + [Poly(sym[0]**N), Poly(sym[1]**N)]
    G_loc = groebner(f_loc, sym, order='grevlex')
    groebner_m_test(4, G_loc, sym)
    Bounds_loc = [g.LM(order='grevlex').exponents for g in G_loc]
    Q_loc = twovar_list(Bounds_loc)
    Q_Mon = [Monomial(exp, sym) for exp in Q_loc]
    #print(Q_loc)
    #E reduced mod G_loc
    E_red = G_loc.reduce(E.as_expr())[1]
    #print("E reduced, ", E_red)
    #Take the first non-zero monomial term in E_red to be the basis
    #element that phi evaluates to be non-zero
    E_mon = E_red.terms()[0]
    #print("dim of Q_0(f)", len(Q_Mon))
    AA = []
    for a in Q_Mon:
        Row = []
        for b in Q_Mon:
            c = a*b
            q = G_loc.reduce(c.as_expr())[1]
            #print("c", c, 'reduced', q)
            #print(div(q, E_red)[0])
            #print(phi(c, E_mon, G))
            Row.append(phi(c, E_mon, G_loc))
        AA.append(Row)
    AA = Matrix(AA)
    #print(AA)
    return similarity_diagonalize(AA)

def three_dim_EKL(f, sym):
    """Input is a polynomial function f, which in this case is a list of
    two polynomials in two variables. Sym is a list of the symbols
    used as variables in the polynomial.

    The output is the diagonal entries in the EKL class
    """
    E = compute_E(f, sym)
    G = groebner(f, sym, order='grevlex')
    if not G.is_zero_dimensional:
        raise TypeError("The polynomial chosen does not yield a 0-dimensional variety")
    #Bounds tells us the outer limits of the basis for the quotient 
    Bounds = [g.LM(order='grevlex').exponents for g in G]
    print(Bounds)
    Quotient = threevar_list(Bounds)
    print(Quotient)
    N = len(Quotient)
    print(N, G)
    #N is the dimension of the quotient ring as a k-vector space. This
    #enables us to sort out the localization at 0 by further modding
    #out.
    f_loc = f + [Poly(sym[0]**N), Poly(sym[1]**N), Poly(sym[2]**N)]
    G_loc = groebner(f_loc, sym, order='grevlex')
    Bounds_loc = [g.LM(order='grevlex').exponents for g in G_loc]
    Q_loc = threevar_list(Bounds_loc)
    print(len(Q_loc), G_loc)
    Q_Mon = [Monomial(exp, sym) for exp in Q_loc]
    #E reduced mod G_loc
    E_red = G_loc.reduce(E.as_expr())[1]
    print("E reduced, ", E_red)
    #Take the first non-zero monomial term in E_red to be the basis
    #element that phi evaluates to be non-zero
    E_mon = E_red.terms()[0]
    #print("dim of Q_0(f)", len(Q_Mon))
    AA = []
    for a in Q_Mon:
        Row = []
        for b in Q_Mon:
            c = a*b
            q = G_loc.reduce(c.as_expr())[1]
            #print("c", c, 'reduced', q)
            #print(div(q, E_red)[0])
            #print(phi(c, E_mon, G))
            Row.append(phi(c, E_mon, G_loc))
        AA.append(Row)
    AA = Matrix(AA)
    #print(AA)
    return similarity_diagonalize(AA)


def rand_poly(deg_max, deg_min, n, sym):
    """gives a random poly of degree deg with integer coefficients between
    -n and n, using the variables of sym.

    """
    m = len(sym)
    L = itermonomials(sym, deg_max, 0)
    out = Poly('0', sym)
    for l in L:
        test = Poly('0', sym)
        if (test + l).total_degree() > deg_min - 1:
            out += randint(-n,n) * l
    return out



def signature(diag):
    out = 0
    for x in diag:
        if x>0:
            out +=1
        elif x<0:
            out -= 1
        else:
            print("diag for signature error", diag)
            raise TypeError
    return out
