import random
import pickle
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit
from qiskit.opflow import X, Z, I, Y
# Use Aer's qasm_simulator
from tqdm import tqdm
from qiskit.opflow import PrimitiveOp

h = np.loadtxt('matrix_Hamiltonian_Bo7',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
h_ex = np.zeros([32,32],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        h_ex[i,j] = h[i,j]

H = PrimitiveOp(h_ex)

def Pauli(i):
    if i == 0:
        return np.array([[1,0],[0,1]])
    if i == 1:
        return np.array([[0,1],[1,0]])
    if i == 2:
        return np.array([[0,-1j],[1j,0]])
    if i == 3:
        return np.array([[1,0],[0,-1]])

def state(x):
    y = np.loadtxt('time_{}_Bo7_new'.format(x),dtype=complex)
    return y


#H_matrix = np.loadtxt('hamiltonian_18F',dtype=complex)
simulator = QasmSimulator()
def evolution(number,Hamiltonian,time,index):
    #for i in tqdm(range(100)):
    circuit = QuantumCircuit(number, number)
    circuit.initialize(state(index))
    list_gate = []
    #list_gate_y = []
    gate_clifford = qiskit.quantum_info.random_clifford(num_qubits=number)
    #print(gate_clifford)
    list_gate.append(gate_clifford.to_operator())
    #print(gate_clifford)
    gate_clifford = gate_clifford.to_instruction()
    circuit.append(gate_clifford,qargs=list(range(number)))
    #print(circuit)
    circuit.measure(list(range(number)), list(range(number)))
    compiled_circuit = transpile(circuit, simulator)
    #for i in tqdm(range(0,100)):
    job = simulator.run(compiled_circuit, shots=1)
    result = job.result()
    counts = result.get_counts(compiled_circuit)
    #print(counts)
    y = list(counts.keys())[list(counts.values()).index(1)]
    state_list = []
    circuit0 = QuantumCircuit(number, number)
    x_clifford = list_gate[0]
    circuit0.append(x_clifford,qargs=list(range(number)))
    z = qiskit.quantum_info.Statevector.from_label(y)
    #print(z)
    #print(circuit0.draw())
    Gate = qiskit.quantum_info.Operator.from_circuit(circuit0)
    Gate = Gate.adjoint()

    z = z.evolve(Gate)

    z_list = []
    for k in range(0, 2**number):
        z_list.append(z[k])
    #print(z_list)
    z_array = np.array([z_list])
    z_matrix = z_array.T.conjugate().dot(z_array)

    z_matrix = z_matrix * (2**number + 1)

    state_list.append(z_matrix)


    return state_list
#print(evolution(6,H2_op,time=1))

#print(evolution(6,H,time=1).draw())


def identity(number):
    matrix = 1
    for i in range(0,number):
        matrix = np.kron(matrix,Pauli(0))
    return matrix
#M_rocover = 3 * Pauli(1) - Pauli(0)

identity10 = identity(5)

def shadow(number,Hamiltonian,time,index):
    state_list = evolution(number,Hamiltonian,time=time,index = index)
    #print(state_list)
    state_list = state_list - identity10
    #print(state_list)

    #print(np.trace(state_matrix))
    return state_list[0]

state_shadow = 0
# for i in tqdm(range(1,21)):
#     print(state(i))
#shadow(10,H,1,index = 1)
for j in tqdm(range(0,100,1)):
    for i in range(0,100):
        state_shadow = (state_shadow + shadow(5,H,1,index = 4)/100)
    if j % 10 == 0:
        np.savetxt('shadow_{}_1h_Bo7_4_94_new_bar_8'.format(j),state_shadow)
#
# number = 10
# simulator = QasmSimulator()
# circuit = QuantumCircuit(number, number)
# circuit.initialize(state(1))
# print(circuit)

#print(state,state+shadow(2,H2_op,0))
# for i in tqdm(range(0,100000)):
#     #pool.apply_async(func=shadow(), args=(i, i + 1))
#     state = (state + shadow(10,H,1)/100000)
    #state0 = state0 + shadow(6,H,0)/1000

