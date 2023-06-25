import numpy as np
from qiskit.opflow import X, Z, I, Y
from tqdm import tqdm


# def Mr():
#     M = np.loadtxt('matrix_Hamiltonian_test')
#     a,b = np.linalg.eig(M)
#     #print(a,b)
#     a = np.delete(a,[-1,-2])
#     #del a[-1]
#     a = np.diag(a)
#     b = np.delete(b,[-1,-2],axis=0)
#     b = np.delete(b,[-1,-2],axis=1)
#     #M_r = (np.linalg.inv(b).dot(a)).dot(b)
#     M_r = (b.dot(a)).dot(np.linalg.inv(b))
#     #a,b = np.linalg.eig(M_r)
#     return M_r


def Mr():
    M = np.loadtxt('test')
    return M

def Pauli(i):
    if i == 0:
        return np.array([[1,0],[0,1]])
    if i == 1:
        return np.array([[0,1],[1,0]])
    if i == 2:
        return np.array([[0,-1j],[1j,0]])
    if i == 3:
        return np.array([[1,0],[0,-1]])


def pauli_index(i):
    if i==0:
        return I
    if i==1:
        return X
    if i==2:
        return Y
    if i==3:
        return Z


M_r = Mr()
# a,b = np.linalg.eig(M_r)
# print(a)

# def Hamiltonian_single(a,b,c,d,e,f):
#     P = np.kron(np.kron(np.kron(np.kron(np.kron(Pauli(a),Pauli(b)),Pauli(c)),Pauli(d)), Pauli(e)),Pauli(f))
#     return  np.trace(M_r.dot(P))/64
#
# #print((0.9 * I^I^I^I^I^Z) + (1.0 * I^I^I^I^I^X))
#
# #print(Hamiltonian_single(1,1,1,1,1,1))
#
# def Hamiltonian():
#     H = 0
#     for a in range(0,1):
#         for b in range(0,1):
#             for c in range(0,1):
#                 for d in range(0,1):
#                     for e in range(0,4):
#                         for f in range(0,4):
#                             H =  H + (Hamiltonian_single(a,b,c,d,e,f)  * pauli_index(a)^pauli_index(b)^pauli_index(c)^pauli_index(d)^pauli_index(e)^pauli_index(f))
#     return H

#print(Hamiltonian())
# H = Hamiltonian()
# a,b = np.linalg.eig(H)
# print(min(a))

def Hamiltonian_single(list):
    P = 1
    for i in list:
        P = np.kron(P,Pauli(i))
    return   np.trace(M_r.dot(P))/64

#print(Hamiltonian_single(0,0,0,0,0,2))


def Hamiltonian():
    H = 0
    for a in range(0,4):
        for b in tqdm(range(0,4)):
            for c in range(0,4):
                for d in range(0,4):
                    for e in tqdm(range(0,4)):
                        for f in range(0,4):
                            H = H + ((Hamiltonian_single([a,b,c,d,e,f])) * pauli_index(a)^pauli_index(b)^pauli_index(c) ^ pauli_index(d) ^ pauli_index(e) ^ pauli_index(f))#^pauli_index(g)^pauli_index(h)^pauli_index(i)^pauli_index(j))
    return H
#np.savetxt('H_pauli_19F',Hamiltonian())



#print(Hamiltonian())

# H = Hamiltonian()
# a,b = np.linalg.eig(H)
# print(min(a))
