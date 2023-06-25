import math

from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit
from qiskit.quantum_info import Statevector
from qiskit.opflow import X, Z, I, Y
# Use Aer's qasm_simulator
from tqdm import tqdm
import matrix
import multiprocessing
from scipy.optimize import fsolve
#import subspace_calculate
import numpy as np
import scipy
import matplotlib.pyplot as plt

#H_matrix = np.loadtxt('hamiltonian_18F',dtype=complex)
h = np.loadtxt('test',dtype=complex)



H_matrix = h

#a,b = np.linalg.eig(H_matrix)
#print(min(a))

# state1 = np.loadtxt('time_1',dtype=complex)
#
# value = state1.dot(H_matrix.dot(state1.T.conjugate()))
#
# print(value)
#
# state2 = np.loadtxt('time_2',dtype=complex)
#
# value = state2.dot(H_matrix.dot(state2.T.conjugate()))
#
# print(value)

def state(x):
    y = np.loadtxt('t={}_state'.format(x),dtype=complex)
    #print(y.shape)
    return y
#z = Statevector(state(3))
#z = qiskit.quantum_info.Statevector.from_label(state(3))
#print(z)

scale = 3
N_matrix = np.zeros((scale,scale),dtype=complex)
M_matrix = np.zeros((scale,scale),dtype=complex)

#print(state(1))

def state_vector0(i):
    # np.random.seed(i)
    # x = np.random.rand(1)
    state_vector = state(i)
    state_vector_list = []
    # print(state_vector)
    for k in range(0, len(state_vector)):
        state_vector_list.append(state_vector[k])
    state_vector = np.array([state_vector_list])
    return state_vector #* (math.e ** (1j * x))


x_value=[]
y_value = []
for n in range(1000,100000,1000):
    for i in range(1, scale + 1):
        for j in tqdm(range(1, scale + 1)):
            # a= H_evolution(10,H,time1=i,time2=j)
            N_matrix[i - 1, j - 1] = state_vector0(i).conjugate().dot(state_vector0(j).T)
            M_matrix[i - 1, j - 1] = state_vector0(i).conjugate().dot(H_matrix.dot(state_vector0(j).T))
    x_value.append(n)
    stadard = 1/np.sqrt(n)
    M_matrix = M_matrix + np.random.normal(0,stadard,[scale,scale]) # + 1j * np.random.normal(0,stadard,[scale,scale])
    #print(M_matrix)
    N_matrix = N_matrix + np.random.normal(0,stadard,[scale,scale]) #+ 1j * np.random.normal(0,stadard,[scale,scale])
    a,b = scipy.linalg.eig(M_matrix,N_matrix)
    value = -abs(min(a))
    a, b = np.linalg.eig(h)
    value_ideal = min(a)
    print(value_ideal,value)
    #y_value.append(stadard)
    y_value.append(1/abs(value_ideal - value))
#np.savetxt('H_19F_matrix_ideal_vector',M_matrix)
#np.savetxt('N_19F_matrix_ideal_vector',N_matrix)

plt.scatter(x_value,y_value)
plt.ylim(0,50)
plt.show()
#
# value = state(3).dot(H_matrix.dot(state(6).T.conjugate()))
#
# print(value)
#
# value = state(3).dot(state(3).T.conjugate())
#
# print(value)