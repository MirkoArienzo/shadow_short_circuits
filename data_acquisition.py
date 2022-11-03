import numpy as np
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import argparse

import symplectic as sy
import shadows as sd
from time import time

parser = argparse.ArgumentParser(description="Pauli estimation simulation")
parser.add_argument("-n", "--nqubits", type=int, default=4, help="Number of qubits")
parser.add_argument("-m", "--samples", type=int, default=10000, help="Number of samples")
parser.add_argument("--pauli", type=str, choices=["X", "Y", "Z"], default="Z", help="Pauli operator")
parser.add_argument("--supp", choices=["half", "full", "threshold"], help="Support of Pauli string")
parser.add_argument("--outfile", type=str, help="Output file")
parser.add_argument("--ensemble", type=str, choices=["BW","LC"], default="BW", help="Ensemble of unitaries")

args = parser.parse_args()
n = args.nqubits
assert n%2==0, "Number of qubits must be even"
N = 2 * n # Phase space dimension
d = 2 ** n # Global dimension
m = args.samples
pauli = args.pauli
supp = args.supp
out = args.outfile
ensemble = args.ensemble


# Generators input state |0>
generators = []
for i in range(n):
    gen = np.zeros(N, dtype = np.uint8)
    gen[2*i] = 1
    generators.append(gen)
    
# Pauli observable described as bitstring in \FF_2^{2n}
if supp=="half":
    v = np.load("./"+str(n)+"_qubits/half_supp_"+str(pauli)+".npy")
if supp=="full":
    v = np.load("./"+str(n)+"_qubits/full_supp_"+str(pauli)+".npy")
if supp=="threshold":
    v = np.load("./"+str(n)+"_qubits/treshold_supp_"+str(pauli)+".npy")
frameop = sd.frameop(v, ensemble)


# print("v = ", v)
# if np.any(sy.lag_x_component_product_coord(v)):
    # print("<0|W(v)|0> = 0")
# else:
    # print("<0|W(v)|0> = 1")
    
    
# print("m = " + str(m))   
outcomes = np.zeros(m) 
outcome = 0
t1 = time()
for j in range(m):
    g = sd.sample_symplectic_part(n, ensemble)
    v2 = sd.update_vector(v, g, ensemble)
    if not np.any(sy.lag_x_component_product_coord(v2)):
        if not np.any(sy.lag_x_component_product_coord(v)):
            # Case v_z = 0, meaning the estimated outcome is 1. We know there cannot be -1 outcomes
            # Therefore, outcome in this case is +1, since we already checked (gv)_x = 0
            outcome += 1
            outcomes[j] = 1
            continue
        a = sd.sample_weyl_part(n, ensemble)
        y = sd.sample_vector(generators, g, a, ensemble)
        sy_form = sd.linear_phase(v, g, a, ensemble)
        dot_prod = sy.dot_product(sy.lag_z_component_product_coord(v2), y)
        alpha = sd.alpha(v, g, ensemble)
        phase = (alpha
                  + sy_form
                  + dot_prod
                  ) % 2
        if phase == 0:
            outcome += 1
            outcomes[j] = 1
        else:
            outcome -= 1
            outcomes[j] = -1
t2 = time()
avg = outcome / (sd.frameop(v, ensemble) * m)
# print("estimated outcome = " + str(avg))
# print("time = ", t2-t1)


path = './'+str(n)+'_qubits/data/'

# Check whether the specified path exists or not
exists = os.path.exists(path)
if not exists:
    # Create a new directory because it does not exist
    os.makedirs(path)

np.save(path+ensemble+'_'+supp+'_'+pauli+'_'+out+'.npy', outcomes)




