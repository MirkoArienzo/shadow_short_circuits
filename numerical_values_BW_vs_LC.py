import numpy as np
import matplotlib.pyplot as plt
from shadows import frameop_periodic, frameop_open 
from matplotlib import rc
import os

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


num_qubs = [4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80]
var_periodic = []
var_open = []
var_LC = []
var_LC_half = []


for n in num_qubs:
    
    # var_periodic = 1 / frameop_periodic(n)
    var_open = 1 / frameop_open(n)
    var_LC = 3**(n)
    var_LC_half = 3**(n/2)
    threshold = 0
    var_LC_threshold = 0
    
    for k in range(n//2, n):
        var_LC_threshold = 3**k
        if var_open <= var_LC_threshold:
            threshold = k
            break
    print("num qubits = " + str(n))
    print("var_BW open = " + str(var_open))
    print("var_LC = " + str(var_LC))
    print("var_LC_half = " + str(var_LC_half))
    print("threshold reached at " + str(k) + " qubits in the support")
    print("var_LC_threshold = " + str(var_LC_threshold))
    print("\n")
    
    

example_open = 1/frameop_open(6) * 1/frameop_open(8) * 1/frameop_open(4)
example_LC = 3**(3) * 3**(6) * 3**(2)
example_open_shifted = 1/frameop_open(6) * 1/frameop_open(6) * 1/frameop_open(6)
print("open = " + str(example_open))
print("LC = " + str(example_LC))
print("shifted = " + str(example_open_shifted))