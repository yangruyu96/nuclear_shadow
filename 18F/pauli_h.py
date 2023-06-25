import numpy as np
from tqdm import tqdm

h = np.loadtxt('hamiltonian_18F',dtype=complex)
h_ex = np.zeros([256,256],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        h_ex[i,j] = h[i,j]

def Pauli(i):
    if i == 0:
        return np.array([[1,0],[0,1]])
    if i == 1:
        return np.array([[0,1],[1,0]])
    if i == 2:
        return np.array([[0,-1j],[1j,0]])
    if i == 3:
        return np.array([[1,0],[0,-1]])

def multi_trace(list):
    A = list[0]
    for i in range(1,len(list)):
        A = A.dot(list[i])
    return np.trace(A)

def Pauli_tensor(x):
    l = len(x)
    p = 1
    for i in range(0,l):
        p = np.kron(p,Pauli(x[i]))
    return p

H_pauli = []
term = 0

for a in tqdm(range(0, 4)):
    for b in range(0, 4):
        for c in range(0, 4):
            for d in range(0, 4):
                for e in range(0, 4):
                    for f in range(0, 4):
                        for g in range(0, 4):
                            for h in range(0, 4):
                                H_pauli.append(multi_trace([h_ex,Pauli_tensor([a,b,c,d,e,f,g,h])]))
                                term = term + 1

np.savetxt('hamiltonian_18F_pauli',H_pauli)
