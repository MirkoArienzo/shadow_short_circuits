"""
    Plot of p(\sigma_{BW}^2 \leq \sigma_{LC}^2)
"""

import numpy as np
import shadows as sd
import matplotlib.pyplot as plt
from matplotlib import rc
import os

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


num_test = 2**16
num_qubs = np.arange(10, 62, 2)
list_ratios = []
avg_differences = []
avg_ratios = []
for n in num_qubs:
    N = 2 * n
    d = 2 ** n
    num_success = 0
    for i in range(num_test):
        v = np.random.randint(0, 2, N, dtype = np.uint8)
        k = 0
        for i in range(0, N, 2):
            if v[i] != 0 or v[i+1] != 0:
                k += 1
        bw_var = 1 / sd.frameop(v, 'BW')
        singlequbs_var = 3**k
        difference = singlequbs_var - bw_var
        ratio = bw_var / singlequbs_var
        if bw_var <= singlequbs_var:
            num_success += 1
    avg_differences.append(difference/num_test)
    avg_ratios.append(ratio**(1/num_test))
    ratio = num_success/num_test
    list_ratios.append(ratio)


plt.figure()
plt.xlabel("$n$", fontsize=16)
plt.ylabel(r"$p(\sigma_{BW}^2 \leq \sigma_{LC}^2)$", fontsize=16)
plt.plot(num_qubs, list_ratios, 'ko')

path = './figures/'
# Check whether the specified path exists or not
exists = os.path.exists('./figures/')
if not exists:
    # Create a new directory because it does not exist
    os.makedirs(path)
plt.savefig('figures/BW_vs_LC_probabilistic.png', dpi = 1200)
plt.show()
