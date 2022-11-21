# Shadow estimation with brickwork circuits 
## [Reference](https://arxiv.org/abs/2211.09835)

Numerical implementation of the Pauli estimation task using different ensembles of unitary operators.

## Scripts:

- In <ins>recursion_open.nb</ins> and <ins>recursion_periodic.nb</ins> the solutions of the systems of recurrence relations are calculated.

Here we mainly compare brickwork with two layers (BW) and local Cliffords (LCs) ensembles, the global Clifford ensembles is implemented for testing only.
The code is self-contained, using standard Python packages only for the simulation scripts.

- `pauli_generator.py': generate three fixed bitstrings given the number of qubits, representing the chosen Pauli observable. The observable is a tensor product of one between X, Y, Z operators and identities.
 The three bitstrings are saved as `.npy' in the `n_qubits' subfolder, where n is the number of qubits, and loaded when needed.
- `data_acquisition.py': collect classical snapshots of an observable of a given type (X,Y,Z), for the given ensemble (LC,BW), for a fixed number of qubits and the chosen bitstring from the ones generated with pauli_generator.py.
- `avg_estimator_comparison.py': generate plots as in [Figure 6 of the paper](https://arxiv.org/abs/2211.09835). In particular, one needs to collect at least 55.000.000 samples to reproduce the figures for each bitstring and ensemble, and save it in the `plot_comparison' subfolder.
- `BW_vs_LC.py': reproduce [Figure 4a of the paper](https://arxiv.org/abs/2211.09835) and save it in the `figures' subfolder.
- `BW_full_supp_vs_LC.py': reproduce [Figure 4a of the paper](https://arxiv.org/abs/2211.09835) and save it in the `figures' subfolder.
- `BW_vs_LC_random.py': reproduce [Figure 5 of the paper](https://arxiv.org/abs/2211.09835) and save it in the `figures' subfolder.
