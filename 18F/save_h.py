from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit
from qiskit.opflow import X, Z, I, Y, OperatorBase, PauliSumOp
# Use Aer's qasm_simulator
from tqdm import tqdm
import pickle
import matrix
import multiprocessing
from scipy.optimize import fsolve
#import subspace_calculate
import numpy as np
import scipy

H_pauli = np.loadtxt('hamiltonian_18F_pauli',dtype=complex)
H_pauli = H_pauli.real
print(H_pauli)
term=0
H = 0
def pauli_index(i):
    if i==0:
        return I
    if i==1:
        return X
    if i==2:
        return Y
    if i==3:
        return Z
for a in tqdm(range(0, 4)):
    for b in range(0, 4):
        for c in range(0, 4):
            for d in range(0, 4):
                for e in range(0, 4):
                    for f in range(0, 4):
                        for g in range(0, 4):
                            for h in range(0, 4):
                                H = H + (H_pauli[term] * pauli_index(a)^pauli_index(b)^pauli_index(c) ^ pauli_index(d) ^ pauli_index(e) ^ pauli_index(f)^pauli_index(g)^pauli_index(h))
                                term = term + 1


# H = np.array([H])
# print(H)
# np.savetxt('Ham',H,fmt = '%s')
#H = OperatorBase(H)
#H = np.loadtxt('Ham',dtype=OperatorBase)
#print(type(H))
#print(H.primitive)
#H0 = PauliSumOp(H.primitive)
#print(type(H0))

output=open('ham_89.pkl','wb')
pickle.dump(H,output)
output.close()

#H_read = PauliSumOp()
# with open('ham.pkl','rb') as file:
#     H_read = pickle.load(file.read())
f = open('ham_89.pkl','rb')
H_read = pickle.load(f)
#H_read.show
print(H_read)
print(type(H_read))