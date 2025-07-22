import numpy as np
from qiskit.circuit import QuantumCircuit, ParameterVector

def makewave0(wavefunction, dele, name):
    n = wavefunction.num_qubits
    param = ParameterVector(name, int(2*n))
    t = 0
    for depth in range(1):
        for i in range(n):
            wavefunction.ry(param[t], i)
            t += 1
        for j in range(int(n/2)):
            wavefunction.cnot(2*j+1, 2*j)
        for i in range(n):
            if i > 0 and i < n-1:
                wavefunction.ry(param[t], i)
                t += 1
        for j in range(int(n/2)):
            if j > 0 and j < int(n/2):
                wavefunction.cnot(2*j, 2*j-1)
    return wavefunction

def makeinitial(wavefunction, param):
    n = wavefunction.num_qubits
    t = 0
    for depth in range(1):
        for i in range(n):
            wavefunction.ry(param[t], i)
            t += 1
        for j in range(int(n/2)):
            wavefunction.cnot(2*j+1, 2*j)
        for i in range(n):
            if i > 0 and i < n-1:
                wavefunction.ry(param[t], i)
                t += 1
        for j in range(int(n/2)):
            if j > 0 and j < int(n/2):
                wavefunction.cnot(2*j, 2*j-1)
    return wavefunction

def makewave1(wavefunction, dele, name):
    n = wavefunction.num_qubits
    param = ParameterVector(name, int(4*n+1))
    t = 0
    for depth in range(1):
        wavefunction.barrier()
        for i in range(n):
            wavefunction.ry(param[t], i)
            t += 1
        for j in range(n-1):
            wavefunction.cnot(j+1, j)
        wavefunction.cnot(0, n-1)
        wavefunction.barrier()
        for i in range(n):
            wavefunction.ry(param[t], i)
            t += 1
        for j in range(n-1):
            wavefunction.cnot(j, j+1)
        wavefunction.cnot(n-1, 0)
    return wavefunction 

# for VITE
def L(params,wavefunction):
    a={}
    t=0
    for i in wavefunction.parameters:
        a[i]=params[t]
        t+=1
        
    qc = wavefunction.assign_parameters(a)
    qc.save_statevector()
    qc.measure_all()
    circ_noise = transpile(qc, sim_noise)
    noise_result = sim_noise.run(circ_noise, shots=1).result()
    u = noise_result.get_statevector()
    return u.conj().dot(Hp.dot(u)).real
def dtheta(params,wavefunction,H):
    N=wavefunction.num_parameters
    A=np.zeros([N,N],dtype=np.complex128)
    C=np.zeros(N,dtype=np.complex128)
    phi=Lv(params,wavefunction)
    dpdt=[]
    cp=1/2
    a=np.pi/2
    for i in range(len(params)):
        ptmp1=params.copy()
        ptmp2=params.copy()
        ptmp1[i]+=a
        ptmp2[i]-=a    
        dp=cp*(Lv(ptmp1,wavefunction)-Lv(ptmp2,wavefunction))
        dpdt.append(dp)
    for i in range(len(params)):
        for j in range(len(params)):
            A[i,j]=(dpdt[i].conj().dot(dpdt[j])).real+dpdt[i].conj().dot(phi)*dpdt[j].conj().dot(phi)
    for i in range(len(params)):
        shape=np.size(dpdt[i])
        # phi=Lv(params,wavefunction)
        C[i]=(dpdt[i].conj().reshape(1,shape).dot((H.dot(phi)).reshape(shape,1))).real
    dx=np.linalg.pinv(A.real).dot(-C)
    return dx.real

def Lv(params,wavefunction):
    a={}
    t=0
    for i in wavefunction.parameters:
        a[i]=params[t]
        t+=1
    qc = wavefunction.assign_parameters(a)
    qc.save_statevector()
    qc.measure_all()
    circ_noise = transpile(qc, sim_noise)
    noise_result = sim_noise.run(circ_noise,shots=1).result()
    u=noise_result.get_statevector()
    return np.array(u)