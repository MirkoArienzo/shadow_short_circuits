"""
    Symplectic tools for shadow estimation, see also  https://kups.ub.uni-koeln.de/50465/
    for theoretical discussion
"""

import numpy as np
from itertools import chain

## return [a, b] = `sum_{i=1}^n a[i]b[i+n] - a[i+n]b[i]
def symplectic_product(a, b):
    n = int(len(a) / 2)
    assert n%2==0, "Vectors must be even dimensional"
    assert len(a)==len(b), "Dimension mismatch"
    a1 = np.array([a[2*i] for i in range(n)])    
    a2 = np.array([a[2*i+1] for i in range(n)])
    b1 = np.array([b[2*i] for i in range(n)])    
    b2 = np.array([b[2*i+1] for i in range(n)])
    return (np.dot(a1, b2) + np.dot(a2, b1)) % 2

## Computes standard dot product of binary vectors modulo 2
def dot_product(v, w):
    assert len(v)==len(w), "Dimension mismatch"
    res = v*w
    return (res.sum()%2)

## Computes standard dot product of binary vectors modulo 4
def dot_product_mod4(v, w):
    assert len(v)==len(w), "Dimension mismatch"
    res = v*w
    return (res.sum()%4)

## Computes the order of the symplectic group Sp(n,q)
def Sp_order(n, q):
    ret = pow(q, n*n)
    for i in range(1,n+1):
        ret *= pow(q, 2*i) - 1
    return  ret

## Returns symplectic form of two binary vectors in zx coordinates 
def symplectic_form(v, w):
    assert hasattr(v, 'shape') and hasattr(w, 'shape'), "symplectic_form:\nArguments do not have shape attribute"
    assert v.shape == w.shape, "symplectic_form:\nDimension mismatch"
    assert v.shape[0] % 2 == 0, "symplectic_form:\nVectors not even dimensional"
    n = v.shape[0] // 2
    return ( dot_product(v[:n], w[n:]) + dot_product(v[n:], w[:n]) ) %2

## Returns symplectic form of two binary vectors in product coordinates
def symplectic_form_product_coord(v, w):
    assert hasattr(v, 'shape') and hasattr(w, 'shape'), "symplectic_form:\nArguments do not have shape attribute"
    assert v.shape == w.shape, "symplectic_form_product_coord:\nDimension mismatch"
    assert v.shape[0] % 2 == 0, "symplectic_form_product_coord:\nVectors not even dimensional"
    n = v.shape[0] // 2
    return ( np.sum( [v[2*i]*w[2*i+1] for i in range(n)] ) + 
            np.sum( [v[2*i+1]*w[2*i] for i in range(n)] ) ) % 2

## Returns a standard basis vector e_i in direction 1 <= i <= n 
def std_basis_vector(i, n):
    assert i <= n and i > 0, "std_basis_vector()\nindex i exceeds bounds"
    return np.array([0]*(i-1)+[1]+[0]*(n-i))


############################
## BASIC MATRIX FUNCTIONS ##
############################

## Returns the identity matrix of order n
def identity(n):
    return np.identity( n, dtype=np.uint8)

## Returns a n x m matrix with zero entries 
def zero_matrix(n, m):
    return np.zeros( (n,m), dtype=np.uint8 )

## Returns a zero vector of size n
def zero_vector(n):
    return np.zeros(n, dtype=np.uint8)

## Returns direct sum of matrices
def direct_sum(*args):
    nargs = len(args)
    for i in range(nargs):
        assert hasattr(args[i], 'shape'), "direct_sum()\nArgument "+str(i)+" does not have a shape attribute"

    shape = ( np.sum( [args[i].shape[0] for i in range(nargs)] ), np.sum( [args[i].shape[1] for i in range(nargs)] ) )
    ret = zero_matrix( *shape )

    row_low_bound = 0
    col_low_bound = 0
    for i in np.arange(0,nargs):
        row_up_bound = row_low_bound + args[i].shape[0]
        col_up_bound = col_low_bound + args[i].shape[1]
        ret[row_low_bound:row_up_bound,col_low_bound:col_up_bound] = args[i]
        row_low_bound = row_up_bound
        col_low_bound = col_up_bound
    return ret

## Row Echelon Form
# Gaussian elimination of a binary matrix
def ref(A):
    n = A.shape[0]  # rows
    m = A.shape[1]  # columns
    row = 0  # row counter
    for col in range(m):  # Go through all columns
        row_pivot = -1
        for k in range(row, n):  # Find pivot element
            if A[k, col] == 1:
                if k != row:  # Swap rows if necessary
                    A[[k]], A[[row]] = A[[row]], A[[k]]
                row_pivot = k  # Set pivot element
                break  # Found pivot element => exit inner for-loop
        if row_pivot == -1:
            if np.array_equal(A[row, :], np.zeros(m, dtype = np.uint8)):  # If entire row is zero
                for k in range(col+1, m):
                    for p in range(row+1, n):  # Find the row with the most left positioned 1.
                        if A[p, k] == 1:  # Swap rows
                            A[[p]], A[[row]] = A[[row]], A[[p]]
                            row_pivot = p  # Set pivot element
                            break  # Found pivot element => exit inner for-loop
                    if row_pivot != -1:
                        break
        if row_pivot != -1:
            for j in range(row_pivot + 1, n):
                if A[j, col] == 1:
                    A[j, :] = (A[j, :] + A[row, :]) % 2  # Add the rows to cancel the 1 below the col element
            row = row + 1
            if row == n:
                break
    return A

## Reduced row echelon form
def rref(matrix):
    h, w = matrix.shape
    row_perm = list(range(h))
    res = matrix
    i = 0
    j = 0
    while i < h and j < w:
        # look for pivot in column j
        i_piv = -1
        for i_piv_cand in range(i, h):
            if res[row_perm[i_piv_cand], j] & 1:
                i_piv = i_piv_cand
                break
        if i_piv == -1:
            # no pivot found, go to next column
            j += 1
        else:
            # we found the next pivot in row `i_piv`. swap this row with the current row (virtually).
            row_perm[i], row_perm[i_piv] = row_perm[i_piv], row_perm[i]
            # now eliminate the pivot from other rows
            
            # if we want to achieve rref, we have to eliminate in ALL other rows
            elim_rows = chain(range(i), range(i+1, h))
            
            # go through the selected rows ...
            for i2 in elim_rows:
                # ... and perform elimination if needed (i.e. if the row has a `1` in the current column)
                if res[row_perm[i2], j] & 1:
                    res[row_perm[i2], j] = 0
                    res[row_perm[i2], j+1:] += res[row_perm[i], j+1:]
            # go to next row and column
            i += 1
            j += 1
    # apply row permutation to memory
    res[:, :] = res[row_perm]
    # return rank and row-reduced matrix
    return res % 2, i

## Extract basis from a set of vectors using reduced row echelon form
def extract_basis(vectors):
    n = len(vectors)
    d = len(vectors[0])
    aux_mat = np.zeros((d, n), dtype = np.uint8)
    for i in range(n):
        aux_mat[:, i] = vectors[i]
    aux_mat, rank = rref(aux_mat)
    aux_mat = aux_mat % 2
    list_pivot_pos = []
    for i in range(d):
        for j in range(i, n):
            if aux_mat[i, j] == 1:
                # print(j)
                list_pivot_pos.append(j)
                break
    return [vectors[i] % 2 for i in list_pivot_pos]


#################################
## SYMPLECTIC MATRIX FUNCTIONS ##
#################################


## Return symplectic unit J 2n x 2n in zx coordinates
def symplectic_J(n):
    J = np.zeros( (2*n, 2*n) )
    J[:n,n:] = identity(n)
    J[n:,:n] = identity(n)
    return J


## Returns symplectic unit J 2n x 2n in product coordinates
def symplectic_J_product_coord(n):
    J = symplectic_J(1)
    return direct_sum( *( [J]*n ) )

## Check if A is symplectic in product coordinates
def is_symplectic_product_coord(A):
    n = A.shape[0] // 2
    J = symplectic_J_product_coord(n)
    return np.array_equal( A.transpose()@J@A % 2, J )

def is_symplectic(A):
    n = A.shape[0] // 2
    J = np.zeros( (2*n, 2*n) )
    J = symplectic_J(n)
    return np.array_equal( A.transpose()@J@A % 2, J )
   

## Returns callable function representing the symplectic transvection t(x) = x + [x,h]h
def transvection(h):
    return ( lambda x: (x + h * symplectic_form(x,h)) % 2 )

## Returns callable function representing the symplectic transvection t(x) = x + [x,h]h
def transvection_product_coord(h):
    return ( lambda x: (x + h * symplectic_form_product_coord(x,h)) % 2 )


###########################################
## WEYL OPERATORS MANIPULATION FUNCTIONS ##
###########################################


def get_basis(v):
    """
        Obtain the basis vectors of the Lagrangian on which v is defined
        Output is a list of basis vectors
    """
    basis = []
    for j in range(len(v)):
        if v[j] != 0:
            vector = np.zeros(len(v), dtype = np.uint8)
            vector[j] = 1
            basis.append(vector)
    return basis

# def eta_product_coord(v, w):
#     """
#         Evaluate function eta(v, w) := v_z * w_x mod 4
#         can be used for phase of Weyl operator W(v) = i^(v_z \cdot v_x) Z(v_z) X(v_x),
#         where v_z \cdot v_x is mod 4
#     """
#     N = len(v)
#     v_z = np.array([v[i] for i in range(0, N, 2)])
#     w_x = np.array([w[i] for i in range(1, N, 2)])
#     return dot_product_mod4(v_z, w_x)

## Evaluate function gamma(v) := v_z * w_x mod 4
def gamma_product_coord(v):
    assert len(v)%2==0, "Vector must be even dimensional"
    n = len(v) // 2
    vz = np.zeros(n, dtype = np.uint8)
    vx = np.zeros(n, dtype = np.uint8)
    for i in range(n):
        vz[i] = v[2*i]
        vx[i] = v[2*i+1]
    return dot_product_mod4(vz, vx)

## Evaluate phase of W(v + w) = (-1)^\beta(v,w) W(v) W(w)
def beta_product_coord(v, w):
    N = len(v)
    assert N%2==0, "Vectors must be even dimensional"
    assert len(v)==len(w), "Dimension mismatch"
    vz = np.zeros(N // 2, dtype = np.uint8)
    wx = np.zeros(N // 2, dtype = np.uint8)
    for i in range(N // 2):
        vz[i] = v[2*i]
        wx[i] = w[2*i+1]
    return ( gamma_product_coord((v+w) % 2) 
            - gamma_product_coord(v)
            - gamma_product_coord(w)
            + 2 * dot_product(vz, wx) ) % 4

def phase_weyls(vectors):
    """
        Compute the phase that is appearing in the product of many Weyl operators
		W(a_1)W(a_2)...W(a_m) = i^{phi(a_1,...,a_m)}W(a_1+...+a_m)
        This phase is given by
		phi(a_1,...,a_m) = sum_{j=1}^{m-1} phi(a_j, sum_{k=j+1}^m a_k)
        
        vectors = list of vectors a_1,..., a_n in F_2^2n
        n = number of qubits
    """
    N = len(vectors[0])
    res = 0
    for j in range(len(vectors)-1):
        temp = np.zeros(N)
        for k in range(j+1, len(vectors)):
            temp = (temp + vectors[k])
        res += beta_product_coord(vectors[j], temp)
    return res % 4

## Check if v is in Lag(x); it is sufficient to check if v commutes with a basis of Lag(x)
def is_in_lagrangian(v, x):
    basis = get_basis(x)
    for basis_vec in basis:
        if symplectic_form_product_coord(v, basis_vec) != 0:
            return False
    return True

def lag_x_component_product_coord(v):
    N = len(v)
    vec = np.zeros(N // 2, dtype = np.uint8)
    for i in range(N // 2):
        vec[i] = v[2*i+1]
    return vec

def lag_z_component_product_coord(v):
    N = len(v)
    vec = np.zeros(N // 2, dtype = np.uint8)
    for i in range(0, N // 2):
        vec[i] = v[2*i]
    return vec

def linear_comb_mod2(vectors, coefficients):
    out = 0
    for i in range(len(coefficients)):
        out += coefficients[i] * vectors[i]
    return out % 2