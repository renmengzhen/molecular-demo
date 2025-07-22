import numpy as np
import itertools
from qiskit.quantum_info import Pauli

def Heisenberg(N, J_x, J_y, J_z):
    J = [J_x, J_y, J_z]
    pp = 0
    H = np.zeros(shape=(2**N, 2**N))
    for p in ['X', 'Y', 'Z']:
        pauli_p = [p*2] + list(itertools.repeat('I', N-2))
        add_p = [p] + list(itertools.repeat('I', N-2)) + [p]
        comb_p = list(set(itertools.permutations(pauli_p, N-1)))
        comb_p.append(list(add_p))
        H_p = np.zeros(shape=(2**N, 2**N))
        for i in range(len(comb_p)):
            pauli_p = ''
            for j in range(len(comb_p[i])):
                pauli_p = pauli_p + comb_p[i][j]
            H_p = H_p + Pauli(pauli_p).to_matrix()
        H = H + J[pp] * H_p
        pp = pp + 1
    
    pauli_z = ['Z'] + list(itertools.repeat('I', N-1))
    comb_z = list(set(itertools.permutations(pauli_z, N)))   
    H_z = np.zeros(shape=(2**N, 2**N))
    for i in range(len(comb_z)):
        pauli_z = ''
        for j in range(len(comb_z[i])):
            pauli_z = pauli_z + comb_z[i][j]
        H_z = H_z + Pauli(pauli_z).to_matrix()
    
    H_fin = -0.5 * (H + H_z)
    return H_fin 