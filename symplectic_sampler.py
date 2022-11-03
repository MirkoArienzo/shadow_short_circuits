"""
    Sample symplectic matrix in Sp(2n,2) using Koenig-Smolin algorithm
    https://arxiv.org/pdf/1406.2170.pdf
"""
import numpy as np
import symplectic as sy


## Find transvection from binary vector x to binary vector y in product coordinates
def find_transvection(x, y):    
    n = x.shape[0] // 2
    assert x.shape == y.shape, "Dimension mismatch"
    assert x.shape[0] % 2 == 0, "Vectors not even dimensional"
    assert not np.array_equal(x, sy.zero_vector(n)) and not np.array_equal(y, sy.zero_vector(n)), "find_transvection():\nOne of the arguments is a zero vector"
    
    #check trivial cases
    if np.array_equal(x, y):
        return ( lambda v: v, lambda v: v )
    if sy.symplectic_form_product_coord(x, y) != 0:
        return ( lambda v: v, sy.transvection_product_coord( (x + y) % 2) )
    
    # check if there is an i s.t. (x_{2i}, x_{2i+1}), (y_{2i}, y_{2i+1}) != 0
    z = sy.zero_vector(2*n)
    for i in range(0, n):
        if (x[2*i] != 0 or x[2*i+1] != 0) and (y[2*i] != 0 or y[2*i+1] != 0):
            # found the pair
            z[2*i] = (x[2*i] + y[2*i]) % 2
            z[2*i+1] = (x[2*i+1] + y[2*i+1]) % 2
            if (z[2*i] + z[2*i+1]) % 2 == 0:  # they were the same so they added to 00
                z[2*i+1] = 1
                if x[2*i] != x[2*i+1]:
                    z[2*i] = 1
            return ( sy.transvection_product_coord( (y - z) % 2 ), sy.transvection_product_coord( (z - x) % 2 ) )       
    
    # didn’t find a pair
    # look for y==00 and x not at first
    for i in range(0, n):
        ii = 2 * i
        if (x[ii] != 0 or x[ii+1] != 0) and (y[ii] == 0 and y[ii+1] == 0):
            # found the pair
            if x[ii] == x[ii+1]:
                z[ii+1] = 1
            else:
                z[ii+1] = x[ii]
                z[ii] = x[ii+1]
            break
    
    # last case, x==00 and y doesn’t
    for i in range(0, n):
        ii = 2 * i
        if (x[ii] == 0 and x[ii+1] == 0) and (y[ii] != 0 or y[ii+1] != 0):
            # found the pair
            if y[ii] == y[ii+1]:
                z[ii+1] = 1
            else:
                z[ii+1] = y[ii]
                z[ii] = y[ii+1]
            break
    return ( sy.transvection_product_coord( (y - z) % 2), sy.transvection_product_coord( (z - x) % 2) )


## Sample a SymplecticMatrix acting on n qubits in product coordinates
## The function draws a random non zero vector to generate a uniformly distributed random symplectic matrix
def symplectic_sampler(n):
        for N in range(2, 2 * n + 1, 2):
            if N == 2:
                g = sy.identity(2)
            else:
                g = sy.direct_sum(sy.identity(2), g)
            
            # sample non null random vector
            v = np.random.randint(0, 2, size = N, dtype = np.uint8)
            while np.array_equal(v, sy.zero_vector(N)):
                v = np.random.randint(0, 2, size = N, dtype = np.uint8)
            
            e1 = sy.std_basis_vector(1, N)
            t = find_transvection(e1, v)

            #define vector e = e1 + sum_{j=2}^{2n} c_j*e_j, with c_j random bit
            e = sy.std_basis_vector(1, N)
            for j in range(2, N):
                e[j] = np.random.randint(0, 2)
            h0 = t[0](t[1](e))

            c = np.random.randint(0, 2)
            if c == 1:
                v = sy.zero_vector(N)
  
            t = ( sy.transvection_product_coord(v), sy.transvection_product_coord(h0) ) + t

            # make g symplectic
            for j in range(N):
                g[j] = t[0](t[1](t[2](t[3](g[j]))))
        return g

