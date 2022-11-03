"""
    Tools for shadow estimation with BW, local Cliffords, and global Cliffords unitaries
    using the frame theory point of view
"""

import numpy as np
import symplectic as sy
import symplectic_sampler as ss

def frameop(v, ensemble):
    if ensemble == 'global':
        return frameop_global(len(v) // 2)
    if ensemble == 'LC':
        return frameop_LC(v)
    if ensemble == 'BW':
        return frameop_BW(v)

## Frame operator of global Cliffords ensemble, given by 1/(2^n+1)
def frameop_global(n):
    return 1/(2**n+1)

## Frame operator of LC ensemble, given by 1/(3^k), where k is the number of qubits supported by v
def frameop_LC(v):
    n = len(v) // 2
    k = 0
    for i in range(n):
        if v[2*i] != 0 or v[2*i+1] != 0:
            k += 1
    return 1/(3**k)

## Frame operator of BW ensemble in the case of fully supported Pauli operators with periodic boundary conditions
def frameop_periodic(n):
    assert n%2==0, "Number of qubits must be even"
    return  ( (np.sqrt(41)+5)**(n//2) + (-1)**(n//2)
                                  *(np.sqrt(41)-5)**(n//2) ) / (5 * np.sqrt(2))**n

## Frame operator of BW ensemble with open boundary conditions on n qubits, n must be even
def frameop_open(n):
    assert n%2==0, "Number of qubits must be even"
    return 5/(2*np.sqrt(41)) * ( (25-3*np.sqrt(41))*(np.sqrt(41)+5)**(n/2) +
                                 (-1)**(n/2+1)*(25+3*np.sqrt(41))*(np.sqrt(41)-5)**(n/2)
                                 ) / ((5*np.sqrt(2))**n)

## Evaluate frame operator of BW ensemble for arbitrary Pauli string v in FF_2^{2n}
def frameop_BW(v):
    n = len(v) // 2
    assert n%2==0, "Number of qubits must be even"
    v = np.roll(v, 2) #needed to make next step easier, since II layer starts from qubit n.1
    is_linear = False
    for i in range(0, 2*n-1, 4):
        if (v[i] == 0 and v[i+1] == 0) and (v[i+2] == 0 and v[i+3] == 0):
            #roll z
            is_linear = True
            v = np.roll(v, 2*n-i)
            break
    if not is_linear:
        return frameop_periodic(n)
    list_lengths = []
    num_pairs = 0
    for i in range(4, 2*n-1, 4):
        if (v[i] == 0 and v[i+1] == 0) and (v[i+2] == 0 and v[i+3] == 0):
            list_lengths.append(num_pairs)
            num_pairs = 0
        else:
            num_pairs += 1
    if num_pairs != 0:
        list_lengths.append(num_pairs)
    if not list_lengths:
        return 1
    else:
        frameop_linear = 1
        for length in list_lengths:
            if length == n // 2:
                return frameop_periodic(n)
            else:
                frameop_linear *= frameop_open(2*(length + 1))
        return frameop_linear


##### Sample random Clifford unitaries, represented in the phase space representation as (g, a), where
##### g is a symplectic matrix in Sp(2n,2), and a is a bitstring of length 2n, n being the number of qubits

## Sample the symplectic matrix on n qubits according to the given ensemble
def sample_symplectic_part(n, ensemble):
    if ensemble == 'global':
        g = ss.symplectic_sampler(n)
    if ensemble == 'LC':
        g = [ss.symplectic_sampler(1) for i in range(n)]
    if ensemble == 'BW':
        g1 = [ss.symplectic_sampler(2) for i in range(n//2)]
        g2 = [ss.symplectic_sampler(2) for i in range(n//2)]
        g = (g1, g2)
    return g

## Sample the Weyl part of a random Clifford, which is represented by a phase space vector a in FF_2^{2n}
def sample_weyl_part(n, ensemble):
    if ensemble == 'BW':
        a1 = np.random.randint(0, 2, 2*n, dtype = np.uint8)
        a2 = np.random.randint(0, 2, 2*n, dtype = np.uint8)
        return (a1, a2)
    else:
        return np.random.randint(0, 2, 2*n, dtype = np.uint8)

def update_vector(v, g, ensemble):
    if ensemble == 'global':
        return upd_vector_global(v, g)
    if ensemble == 'LC':
        return upd_vector_LC(v, g)
    if ensemble == 'BW':
        return upd_vector_BW(v, g)

def upd_vector_global(v, g):
    return g @ v % 2

def upd_vector_LC(v, g):
    u = np.zeros(len(v), dtype = np.uint8)
    for i in range(len(g)):
        u[2*i:2*i+2] = g[i] @ v[2*i:2*i+2]
    return u % 2

def apply_layer1(vec, g):
    """
        Apply symplectic part of U1 to vec according to BW structure.
        
        input g list of n/2 matrices in Sp(4,2)
    """
    c = 0
    N = len(vec)
    n = len(vec) // 2
    res = np.zeros(N, dtype = np.uint8)
    for j in range(0, n // 2):
        res[4*j:4*j+4] = g[c] @ vec[4*j:4*j+4] % 2
        c += 1
    return res

def apply_layer2(vec, g):
    """
        Apply symplectic part of U2 to vec according to the representation of U2
        
        input g list of n/2 matrices in Sp(4,2).
        Last element of list of sympl is the symplectic matrix acting on (n, 1) qubits
    """
    c = 0
    N = len(vec)
    n = N // 2
    out = np.zeros(N, dtype = np.uint8)
    for j in range(0, n // 2 - 1):
        out[4*j+2:4*j+6] = g[c] @ vec[4*j+2:4*j+6] % 2
        c += 1
    ## Update for n-th and first qubit according the last element
    #  of the list sympl
    temp = np.zeros(4, dtype = np.uint8)
    temp[:2] = vec[N-2:]
    temp[2:] = vec[:2]
    temp = g[n//2 - 1] @ temp % 2
    out[N-2:] = temp[:2]
    out[:2] = temp[2:]
    return out

def upd_vector_BW(v, g):
    ### Apply U_2
    out = apply_layer2(v, g[1])
    ### Apply U_1
    out = apply_layer1(out, g[0])
    return out

def upd_generators(generators, g, ensemble):
    if ensemble == 'global':
        return upd_generators_global(generators, g)
    if ensemble == 'LC':
        return upd_generators_LC(generators, g)
    if ensemble == 'BW':
        return upd_generators_BW(generators, g)

def upd_generators_global(generators, g):
    return [update_vector(gen, g, 'global') for gen in generators]

def upd_generators_LC(generators, g):
    """
        generators is list of n vectors in FF_2^{2n}
        g is list of n 2x2 symplectic matrices
    """
    return [update_vector(gen, g, 'LC') for gen in generators]

def upd_generators_BW(generators, g):
    """
    generator is list of n vectors in FF_2^{2n}
    g = (g1, g2) tuple of 2 lists of n/2 4x4 symplectic matrices
    """
    return [update_vector(gen, g, 'BW') for gen in generators]


def sample_vector(generators, g, a, ensemble):
    upd_gens = upd_generators(generators, g, ensemble)
    proj_vectors = [ sy.lag_x_component_product_coord(gen) for gen in upd_gens]
    proj_basis = sy.extract_basis(proj_vectors)
    coefficients = np.random.randint(0, 2, len(proj_basis), dtype = np.uint8)
    if ensemble == 'BW':
        aff_shift = (sy.lag_x_component_product_coord(a[0]) + sy.lag_x_component_product_coord(a[1])) % 2
    else:
        aff_shift = sy.lag_x_component_product_coord(a)
    return ( sy.linear_comb_mod2(proj_basis, coefficients) + aff_shift ) % 2

def alpha(v, g, ensemble):
    if ensemble == 'global' or ensemble == 'LC':
        return alpha_one_layer(v, g, ensemble)
    if ensemble == 'BW':
        return alpha_BW(v, g[0], g[1])

def alpha_one_layer(v, g, ensemble):
    supp_v = np.flatnonzero(v)
    if len(supp_v) < 2:
        return 0
    i = supp_v[0]
    u = np.copy(v)
    w = np.zeros(len(v), dtype = np.uint8)
    u[i] = 0
    w[i] = 1
    u2 = update_vector(u, g, ensemble)
    w2 = update_vector(w, g, ensemble)
    beta_diff = ( sy.beta_product_coord(u2, w2) - sy.beta_product_coord(u, w) ) // 2
    return (beta_diff + alpha_one_layer(u, g, ensemble)) % 2

def alpha_one_layer_bw(v, g, layer):
    supp_v = np.flatnonzero(v)
    if len(supp_v) < 2:
        return 0
    i = supp_v[0]
    u = np.copy(v)
    w = np.zeros(len(v), dtype = np.uint8)
    u[i] = 0
    w[i] = 1
    if layer == 'first':
        u2 = apply_layer1(u, g) 
        w2 = apply_layer1(w, g)
    else:
        u2 = apply_layer2(u, g)
        w2 = apply_layer2(w, g)
    beta_diff = ( sy.beta_product_coord(u2, w2) - sy.beta_product_coord(u, w) ) // 2
    return (beta_diff + alpha_one_layer_bw(u, g, layer)) % 2

def alpha_BW(v, g1, g2):
    alpha = alpha_one_layer_bw(v, g2, 'second')
    alpha = (alpha + alpha_one_layer_bw(apply_layer2(v, g2), g1, 'first')) % 2
    return alpha

def linear_phase(v, g, a, ensemble):
    if ensemble == 'LC' or ensemble == 'global':
        return sy.symplectic_form_product_coord(update_vector(v, g, ensemble), a)
    else:
        v2 = apply_layer2(v, g[1])
        return (sy.symplectic_form_product_coord(a[1], v) + 
                sy.symplectic_form_product_coord(a[0], v2)
                ) % 2