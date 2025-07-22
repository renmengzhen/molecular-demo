import numpy as np
import itertools
import random
from qiskit.quantum_info import Pauli

def spin(N, seed):
    random.seed(seed)
    J = []
    for pp in range(int(N*(N-1)*1.5)):
        J.append(random.uniform(-1,1))
    t=0
    H_p = np.zeros(shape=(2**N,2**N))
    for p in ['X','Y','Z']:
        pauli_p = [p]+[p]+list(itertools.repeat('I',N-2))
        comb_p = list(set(itertools.permutations(pauli_p,N)))
        comb_p.sort()
        for i in range(len(comb_p)):
            pauli_p = ''
            for j in range(len(comb_p[i])):
                pauli_p = pauli_p+comb_p[i][j]
            H_p = H_p+J[i+t*len(comb_p)]*Pauli(pauli_p).to_matrix()
        t=t+1
    return H_p 