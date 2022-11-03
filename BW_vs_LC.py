import numpy as np
import matplotlib.pyplot as plt
from shadows import frameop_periodic, frameop_open 
from matplotlib import rc
import os

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


num_qubs = np.arange(4, 32, 2)
N = len(num_qubs)
var_periodic = []
var_open = []
var_LC = []
var_LC_half = []

for n in num_qubs:
    n = int(n)
    var_periodic.append(1 / frameop_periodic(n))
    var_open.append(1 / frameop_open(n))
    var_LC.append(3**(n))
    var_LC_half.append(3**(n//2))
          
fig = plt.figure()

ax = fig.add_axes([0.11,0.11,0.75,0.75]) # axis starts at 0.11, 0.11
# ax.set_title("Variance scaling for BW and LC")
ax.set_yscale('log')

default_x_ticks = range(1, len(num_qubs)+1)


# plt.title("Variance scaling")
plt.xlabel("$n$", fontsize=16)
plt.ylabel("$\sigma^2$", fontsize=16)

overlapping = 0.750

plt.plot(num_qubs, var_LC, 'm^', label="LC k = n", alpha = overlapping)
plt.plot(num_qubs, var_LC_half, 'gv', label="LC k = n/2", alpha = overlapping)
plt.plot(num_qubs, var_periodic, 'bo', label="BW periodic", alpha = overlapping)
plt.plot(num_qubs, var_open, 'ro', label="BW open (fully supported)", alpha = overlapping)

plt.legend()

path = './figures/'
# Check whether the specified path exists or not
exists = os.path.exists('./figures/')
if not exists:
    # Create a new directory because it does not exist
    os.makedirs(path)
plt.savefig('./figures/BW_vs_LC.png', dpi = 1200)
plt.show()
