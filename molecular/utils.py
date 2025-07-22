import numpy as np
from qiskit.quantum_info import Pauli, state_fidelity

def commutator(A, B):
    return A.dot(B) - B.dot(A)

def anticommutator(A, B):
    return A.dot(B) + B.dot(A)

def label2Pauli(s):
    xs = []
    zs = []
    label2XZ = {'I': (0, 0), 'X': (1, 0), 'Y': (1, 1), 'Z': (0, 1)}
    for c in s[::-1]:
        x, z = label2XZ[c]
        xs.append(x)
        zs.append(z)
    return Pauli(z=zs, x=xs)

def Tmax(u):
    v = u.copy()
    for i in range(len(v)):
        v[i] = v[i] % (4 * np.pi)
    return v

def normalization(data):
    a = np.size(data, 0)
    b = np.size(data, 1)
    data1 = np.zeros(shape=(a, b))
    data2 = np.zeros(shape=(a, b))
    for i in range(a):
        for j in range(b):
            data1[i, j] = (np.exp(1j * data[i, j])).real
            data2[i, j] = (np.exp(1j * data[i, j])).imag
    return np.hstack((data1, data2))

def sf(v, u, w):
    sf_v = []
    us = np.argsort(u)
    for i in us:
        sf_v.append(state_fidelity(v, w[:, i]))
    return sf_v 