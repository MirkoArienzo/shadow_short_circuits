import numpy as np
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import argparse

import symplectic as sy
import shadows as sd
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


parser = argparse.ArgumentParser(description="errors fit for given pauli observable specified by type and support")
parser.add_argument("-n", "--nqubits", type=int, default=4, help="Number of qubits")
parser.add_argument("--supp", type=str, choices=["half","full","threshold"], help="Support Pauli string")
parser.add_argument("--pauli", type=str, choices=["X","Y","Z"], default="Z", help="Type Pauli observable")


args = parser.parse_args()
n = args.nqubits
assert n%2==0, "Number of qubits must be even"
supp = args.supp
pauli = args.pauli

# Pauli observable described as bitstring in \FF_2^{2n}
if supp=="half":
    v = np.load("./"+str(n)+"_qubits/half_supp_"+str(pauli)+".npy")
if supp=="full":
    v = np.load("./"+str(n)+"_qubits/full_supp_"+str(pauli)+".npy")
if supp=="threshold":
    v = np.load("./"+str(n)+"_qubits/threshold_supp_"+str(pauli)+".npy")
print(v)
if np.any(sy.lag_x_component_product_coord(v)):
    target_outcome = 0
else:
    target_outcome = 1
    
data_path = "./"+str(n)+"_qubits/data/"
data_files = os.listdir(data_path)    
frames = ["BW", "LC"]
colors = ["b", "r"]
color_index = 0
num_samples = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]

fig = plt.figure()

# Creating axes instance
ax = fig.add_axes([0.14,0.14,0.8,0.8])

# Creating plot
overlapping = 0.250
default_x_ticks = range(1, len(num_samples)+1)
plt.xlabel("$m$", fontsize=16)
plt.ylabel(r"$\hat w(v)$", fontsize=16)


for frame in frames:
    if frame=="BW":
        ensemble = "brickwork"
    else:
        ensemble = "local"
    frameop = sd.frameop(v, ensemble)
    var = 1 / frameop
    
    data = np.array([])
    path_to_find = f"{frame}_{supp}_{pauli}"
    for data_file in data_files:
        if path_to_find in data_file:
            data_string = f"./{n}_qubits/data/{data_file}"
            data = np.append(data, np.load(data_string))
    estimators = []
    errors = []
    var_estimators = []
    ind_used_samples = 0
    i = 1
    num_rep = 100
    for mm in num_samples:
        estimator = data[ind_used_samples:ind_used_samples+num_rep*mm].sum() / (frameop*num_rep*mm)
        ind_used_samples += num_rep*mm
        estimators.append(estimator)
        var_estimators.append(np.sqrt(var / (num_rep*mm)))
        
    
    plt.errorbar(num_samples, estimators, yerr = var_estimators, fmt='o', 
                 color=colors[color_index], capsize = 5, alpha = overlapping,
                 label = f"{frame} estimator")
    color_index += 1
    
plt.axhline(y = 1, color = 'k', linestyle = '-', label = "Target outcome")
plt.legend()

save_path = "./plot_comparison/"
exists = os.path.exists(save_path)
if not exists:
    os.makedirs(save_path)
plot_name = f"/{n}_{supp}_{pauli}_estimator_plot.png"
plt.savefig(save_path+plot_name, dpi = 1200)

# plt.show()







