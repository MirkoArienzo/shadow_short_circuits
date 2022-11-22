import numpy as np
import argparse
import os
import shadows as sd


def assign_pauli(pauli):
    v = np.zeros(2, dtype = np.uint8)
    if pauli == "Z":
        v[0] = 1
    else:
        if pauli == "X":
            v[1] = 1
        elif pauli == "Y":
            v[0] = 1
            v[1] = 1
    return v

parser = argparse.ArgumentParser(description="Generate Pauli string of a given type with one qubit supported for each brick")
parser.add_argument("-n", "--nqubits", type=int, default=4, help="Number of qubits (must be even)")
parser.add_argument("--pauli", type=str, default="Z", choices=["X","Y","Z"], help="Local Pauli")

args = parser.parse_args()
n = args.nqubits
assert n%2==0, "Number of qubits must be even"
pauli = args.pauli

N = 2 * n

### As a convention, we choose the first qubit of the brick to be supported.
### This means that for the brick acting on  qubits (n, 1) the n-th qubit is supported
### The string is generated as if it looked at the first layer first, than rotated to make it
### consistent with the second layer of the circuit and the way the frame operator is calculated

## Create subfolder for n qubits if does not exists
folder_path = './'+str(n)+'_qubits'
# Check whether the specified path exists or not
exists = os.path.exists(folder_path)
if not exists:
    # Create a new directory because it does not exist
    os.makedirs(folder_path)

savepath = "./"+str(n)+"_qubits/n.npy"
np.save(savepath, n)

v = np.zeros(N, dtype = np.uint8)

for i in range(N//4):
    v[4*i:4*i+2] = assign_pauli(pauli)

v = np.roll(v, -2)

savepath = "./"+str(n)+"_qubits/half_supp_"+str(pauli)+".npy"
np.save(savepath, v)


### Full support generation
### In this case there is no need of a rotation or of any convention, 
### since every qubit is in the support of v

v = np.zeros(N, dtype = np.uint8)


for i in range(n):
    v[2*i:2*i+2] = assign_pauli(pauli)
    

savepath = "./"+str(n)+"_qubits/full_supp_"+str(pauli)+".npy"
np.save(savepath, v)

### Generation of bitstring with support above threshold
### As a convention, we choose the first qubit of the brick to be supported.
### This means that for the (n, 1) the n-th qubit is supported
### Then, we add more qubits in the support of, until we reach n/2+s+1, where s is the threshold
### where BW starts to outperform LC

v = np.zeros(N, dtype = np.uint8)

for i in range(N//4):
    v[4*i:4*i+2] = assign_pauli(pauli)
    if pauli == "Z":
        v[4*i] = 1
    else:
        if pauli == "X":
            v[4*i+1] = 1
        elif pauli == "Y":
            v[4*i] = 1
            v[4*i+1] = 1

v = np.roll(v, -2)

s = 0
ind = 0
flag = False
for i in range(N//4):
    v[4*i:4*i+2] = assign_pauli(pauli)
    s += 1
    k = 0
    for j in range(0, N, 2):
        if v[j] != 0 or v[j+1] != 0:
            k += 1
    LC_var = 3**k
    v2 = np.roll(v, -2)
    BW_var = 1 / sd.frameop(v, "BW")
    if BW_var <= LC_var:
        ind = i+1
        flag = True
        break

if n!=4 and n!=6 and flag:
    v[4*ind:4*ind+2] = assign_pauli(pauli)
    

savepath = "./"+str(n)+"_qubits/treshold_supp_"+str(pauli)+".npy"
np.save(savepath, v)




