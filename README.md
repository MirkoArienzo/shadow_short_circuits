# Shadow estimation with brickwork circuits 

Numerical implementation of Pauli estimation using short random quantum circuits.

## Table of Contents
- [Overview](#overview)
- [System requirements](#system-requirements)
- [Scripts](#scripts)
- [Quick Start](#quick-start-sample)


# Overview
*shadow\_short\_circuits* is a tiny toolbox to simulate **classical shadows** on *shallow* quantum circuits. 
The goal is to provide a minimal, transparent reference implementation for Ref. https://arxiv.org/abs/2211.09835.

# System requirements
## Hardware requirements
A standard computer with enough RAM to support the in-memory operations is enough to run the simulation.

## Software requirements

### OS requirements
The toolbox should work on every platform.
The package has been tested on the following system:
+ Windows 10 using `python3.9.18`

### Python dependencies
+ numpy>=1.24
+ matplotlib>=3.7

### Mathematica dependencies
The notebooks use only built-in Wolfram language functionality—no add-on packages required
+ Recommended version >= 13.0


# Scripts 
## Python
Compute classical shadows using local Cliffords (LCs) and short brickwork (BW).
Global Clifford unitaries are also implemented for testing.
Pauli operators are considered in their phase space representation.

- <ins>pauli\_generator.py</ins>: generate three fixed bitstrings given the number of qubits. The observable is a tensor product of one between X, Y, Z operators and identities.
 The bitstrings will be stored for further processing.
- <ins>data\_acquisition.py</ins>: collect classical snapshots of an observable of a given type (X,Y,Z), for the given ensemble (LC,BW), for a fixed number of qubits and the chosen bitstring.
- <ins>avg\_estimator\_comparison.py</ins>: generate plots of Figure 6.
- <ins>BW\_vs\_LC.py</ins>: reproduce Figure 4b.
- <ins>BW\_full\_supp\_vs\_LC.py</ins>: reproduce Figure 4b.
- <ins>BW\_vs\_LC\_random.py</ins>: reproduce Figure 5.

## Mathematica
Solve the system of recurrence relations for the eigenvalues of the frame operator.
- <ins>recursion\_open.np</ins>: system of equations for open boundary conditions is solved.
- <ins>recursion\_periodic.np</ins>: system of equations for periodic boundary conditions is solved.

# Quick start sample

+ Generate Pauli strings for 8 qubits
```
python pauli_generator.py --n 8
```
+ Collect 10000 samples using BW unitaries
```
python data_acquisition.py --observable Z --ensemble BW --n 8 --shots 10000
```
+ Plot Fig 4a
```
python BW_vs_LC.py --n 8 --outfile figures/figure4a.png
```
