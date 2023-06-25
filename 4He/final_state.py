from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit
from qiskit.opflow import X, Z, I, Y, OperatorBase, PauliSumOp
# Use Aer's qasm_simulator
from tqdm import tqdm
import pickle
import numpy as np
import scipy
from qiskit.opflow import PrimitiveOp
#f = open('Bo7.pkl','rb')
#H_read = pickle.load(f)
#H_read.show
# print(H_read)
# print(type(H_read))

h = np.loadtxt('matrix_Hamiltonian_He6',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
h_ex = np.zeros([16,16],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        h_ex[i,j] = h[i,j]

H_read = PrimitiveOp(h_ex)

def H_evolution(number,Hamiltonian,time1):
    #istate = np.loadtxt('Li8_initial',dtype=complex)
    #print(istate)
    state_initial = qiskit.quantum_info.Statevector.from_label('0000')
    #print(state_initial)
    simulator = QasmSimulator()
    circuit = QuantumCircuit(number, number)
    #circuit.cx(0, 1)
    #print(list(range(number)))
    #print(circuit.draw())
    #circuit.h(0)
    #circuit.s(1)
    #circuit.cx(0, 1)
    circuit.hamiltonian(Hamiltonian, time=time1, qubits=list(range(number)))
    #circuit1 = QuantumCircuit(number,number)
    #circuit3 = QuantumCircuit(number, number)
    #circuit0.h(0)
    #circuit0.s(1)
    #circuit0.cx(0, 1)
    #circuit1.hamiltonian(Hamiltonian,time=time2,qubits=list(range(number)))
    #circuit3.hamiltonian(Hamiltonian, time=3, qubits=list(range(number)))
    #circuit.h(0)
    #circuit.s(1)
    Gate = qiskit.quantum_info.Operator.from_circuit(circuit)
    #Gate1 = qiskit.quantum_info.Operator.from_circuit(circuit1)
    #Gate3 = qiskit.quantum_info.Operator.from_circuit(circuit3)
    #print(circuit)
    #print(circuit)
    #print(H)
    #print(Gate)
    state_final = state_initial.evolve(Gate)
    #print(state_final)
    #state_final1 = state_initial.evolve(Gate1)
    #state_final3 = state_initial.evolve(Gate3)
    #print(state_final0)
    #print(state_final - state_initial)
    #expectation = state_final.expectation_value(Hamiltonian)
    z = state_final
    z_list = []
    for k in range(0, 2 ** number):
        #print(k)
        z_list.append(z[k])
    # print(z_list)
    z_array = np.array([z_list])
    #z1 = state_final1
    # z1_list = []
    # for k in range(0, 2 ** number):
    #     z1_list.append(z1[k])
    # print(z_list)
    #z1_array = np.array([z1_list])
    #z3 = state_final3
    #z3_list = []
    #for k in range(0, 2 ** number):
    #    z3_list.append(z3[k])
    # print(z_list)
    #z3_array = np.array([z3_list])
    #print(z1_array.dot(z_array.T.conjugate())*z_array.dot(z3_array.T.conjugate())*z3_array.dot(z1_array.T.conjugate()))
    #print(z_array.dot(H_matrix.dot(z3_array.conjugate().T)) * z3_array.dot(z_array.conjugate().T))
    #print(z1_array.dot(z_array.T.conjugate()))
    #print(z)4
    return z_array

a= H_evolution(4,H_read,time1=10 )
np.savetxt('time_10_He6_new',a)
#print(a)
#print(type(a))