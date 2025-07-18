import numpy as np
import shadows as sd
import matplotlib.pyplot as plt
from matplotlib import rc
import os

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)



num_qubits = np.arange(6, 102, 6)
thresholds_periodic = []
thresholds_open = []
num_qubs_bound = np.arange(2, 102, 2)
bound = []

for n in num_qubs_bound:
    bound.append(0.68*(1+1/n))
    
for n in num_qubits:
    N = 2 * n

    v = np.zeros(N, dtype = np.uint8)
    for i in range(N//4):
        v[4*i] = 1
    v = np.roll(v, -2)
    
    count = 0
    for i in range(N//4):
        v[4*i] = 1
        count += 1
        k = 0
        for j in range(0, N, 2):
            if v[j] != 0 or v[j+1] != 0:
                k += 1
        LC_var = 3**k
        v2 = np.roll(v, -2)
        BW_var = 1 / sd.frameop(v, "BW")
        if BW_var <= LC_var:
            thresholds_periodic.append(1/2+count/n)
            break
    

fig = plt.figure()

ax = fig.add_axes([0.11,0.11,0.75,0.75]) # axis starts at 0.1, 0.1

default_x_ticks = range(1, len(num_qubits)+1)
default_y_ticks = range(1, len(thresholds_periodic)+1)
plt.xlabel("$n$", fontsize=16)
plt.ylabel("critical fraction of qubits", fontsize=14)
overlapping = 0.250
plt.plot(num_qubits, thresholds_periodic, 'bo', markersize=3, label="Numerical data")
plt.plot(num_qubs_bound, bound, 'r', label=r'$0.68(1+1/n)$', alpha=0.5)

plt.rcParams['lines.color'] = 'g'
plt.rcParams['lines.linewidth'] = 1.0
# plt.vlines(num_qubits, 0.67, thresholds_periodic, linestyles='dashed')

plt.legend(loc="upper right")

path = './figures/'
# Check whether the specified path exists or not
exists = os.path.exists('./figures/')
if not exists:
    # Create a new directory because it does not exist
    os.makedirs(path)
plt.savefig("./figures/thresholds_BW_LC_normalized.png", dpi=1200)

plt.show()

