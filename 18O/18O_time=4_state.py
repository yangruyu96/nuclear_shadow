import random

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit
from qiskit.opflow import X, Z, I, Y
# Use Aer's qasm_simulator
from tqdm import tqdm
import matrix
import multiprocessing
from scipy.optimize import fsolve

def Pauli(i):
    if i == 0:
        return np.array([[1,0],[0,1]])
    if i == 1:
        return np.array([[0,1],[1,0]])
    if i == 2:
        return np.array([[0,-1j],[1j,0]])
    if i == 3:
        return np.array([[1,0],[0,-1]])

H = matrix.Hamiltonian()
H_matrix = np.loadtxt('test')
# H2_op = (-1.052373245772859 * I ^ I^I) + \
#         (0.39793742484318045 * I ^ Z^I) + \
#         (-0.39793742484318045 * Z ^ I^Z) + \
#         (-0.01128010425623538 * Z ^ Z^Z) + \
#         (0.18093119978423156 * X ^ X^Z)

#H2_op = (I ^ X) #+ (0.18093119978423156 * X ^ X)
# H2_op = I
# H2_op_matrix = Pauli(0)
#
# H2_op_matrix = (-1.052373245772859 * np.kron(np.kron(Pauli(0),Pauli(0)),Pauli(0))) + (0.39793742484318045 *np.kron(np.kron(Pauli(0),Pauli(3)),Pauli(0))) + \
# (-0.39793742484318045 * np.kron(np.kron(Pauli(3),Pauli(0)),Pauli(3))) + (-0.01128010425623538 * np.kron(np.kron(Pauli(3),Pauli(3)),Pauli(3))) + \
# (0.18093119978423156 * np.kron(np.kron(Pauli(1),Pauli(1)),Pauli(3)))
#
# H2_op_matrix = (np.kron(Pauli(2),Pauli(3))) + (0.18093119978423156 * np.kron(Pauli(1),Pauli(3))) + np.kron(Pauli(3),Pauli(0))
# H2_op = (Y^Z) + (0.18093119978423156 * X^Z) + (Z^I)
#print(H2_op_matrix)
#
# H2_op = X^Z
# H2_op_matrix = np.kron(Pauli(1),Pauli(3))
# print(H2_op_matrix)
# H2_op = (-1.052373245772859 * I)+ (0.39793742484318045 * X) + \
# (-0.39793742484318045 * Y) + (-0.01128010425623538 * Z)
#
# H2_op_matrix = -1.052373245772859 * Pauli(0)+ (0.39793742484318045 * Pauli(1)) + \
# (-0.39793742484318045 * Pauli(2)) + (-0.01128010425623538 * Pauli(3))


#H2_op_matrix = np.kron(Pauli(1),Pauli(2))

#print(H_matrix)
def H_evolution(number,Hamiltonian,time):
    state_initial = qiskit.quantum_info.Statevector.from_label('000000' )
    #print(state_initial)
    simulator = QasmSimulator()
    circuit = QuantumCircuit(number, number)
    #circuit.cx(0, 1)
    #print(list(range(number)))
    #print(circuit.draw())
    #circuit.h(0)
    #circuit.s(1)
    #circuit.cx(0, 1)
    circuit.hamiltonian(Hamiltonian, time=time, qubits=list(range(number)))
    circuit0 = QuantumCircuit(number,number)
    #circuit0.h(0)
    #circuit0.s(1)
    #circuit0.cx(0, 1)
    #circuit0.hamiltonian(Hamiltonian,time=1,qubits=list(range(number)))
    #circuit.h(0)
    #circuit.s(1)
    Gate = qiskit.quantum_info.Operator.from_circuit(circuit)
    Gate0 = qiskit.quantum_info.Operator.from_circuit(circuit0)
    #print(circuit)
    #print(circuit)
    #print(H)
    #print(Gate)
    state_final = state_initial.evolve(Gate)
    #print(state_final)
    state_final0 = state_initial.evolve(Gate0)
    #print(state_final0)
    #print(state_final - state_initial)
    expectation = state_final.expectation_value(Hamiltonian)

    return expectation,state_final.inner(state_final0)*state_final.inner(state_final0).conjugate(),state_final
    # simulator = QasmSimulator()
    # circuit = QuantumCircuit(number, number)
    # circuit.hamiltonian(Hamiltonian, time=time, qubits=list(range(number)))
    # #circuit.measure(list(range(number)), list(range(number)))
    # compiled_circuit = transpile(circuit, simulator)
    # #job = simulator.run(compiled_circuit, shots=1000)
    # circuit.save_expectation_value(operator=Hamiltonian,qubits=list(range(0,number)))
    # print(circuit)
    #result = job.result()
    #counts = result.get_counts(compiled_circuit)

#print(H_evolution(6,H,4)) #0.173 0.875 # 0.9676 -1.06365

#print(H_evolution(6,H,1)) #-7.9685

def evolution(number,Hamiltonian,time):
    simulator = QasmSimulator()
    circuit = QuantumCircuit(number, number)
    #circuit.cx(0, 1)
    circuit.hamiltonian(Hamiltonian, time=time, qubits=list(range(number)))
    #circuit.h(0)
    #circuit.x(0)
    #circuit.x(1)
    #circuit.s(1)
    #circuit.cx(0, 1)
    list_gate = []
    list_gate_y = []
    gate_clifford = qiskit.quantum_info.random_clifford(num_qubits=number)
    list_gate.append(gate_clifford.to_operator())
    gate_clifford = gate_clifford.to_instruction()
    circuit.append(gate_clifford,qargs=list(range(number)))
    #print(circuit)
    circuit.measure(list(range(number)), list(range(number)))
    compiled_circuit = transpile(circuit, simulator)
    job = simulator.run(compiled_circuit, shots=1)
    result = job.result()
    counts = result.get_counts(compiled_circuit)
    #print(counts)
    y = list(counts.keys())[list(counts.values()).index(1)]
    #print(y)
    #print(y)
    #print(y[0],y[1])
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
    #print(z)
    z_list = []
    for k in range(0, 2**number):
        z_list.append(z[k])
    #print(z_list)
    z_array = np.array([z_list])
    z_matrix = z_array.T.conjugate().dot(z_array)
    #print(np.trace(z_matrix))
    z_matrix = z_matrix * (2**number + 1)
    #print(z_matrix)
    state_list.append(z_matrix)
    #print(z)
    #print(state_list)
    # print(z_matrix)
    # #z_array = z_array.reshape(number,number)
    # #print(z_array)
    # print(np.trace(z_matrix))
    # print(z.expand())
    #print(z.expectation_value(Hamiltonian))
    #print(counts)
    return state_list
#print(evolution(6,H2_op,time=1))

#print(evolution(6,H,time=1).draw())


def identity(number):
    matrix = 1
    for i in range(0,number):
        matrix = np.kron(matrix,Pauli(0))
    return matrix
#M_rocover = 3 * Pauli(1) - Pauli(0)

def shadow(number,Hamiltonian,time):
    state_list = evolution(number,Hamiltonian,time=time)
    #print(state_list)
    state_list = state_list - identity(6)
    #print(state_list)

    #print(np.trace(state_matrix))
    return state_list[0]
    #print(state_list)
    #print(np.trace(state_matrix))
    #print(state_matrix - state_matrix.T.conjugate())
    #for i in range(0,number):

#pool = multiprocessing.Pool(8)
# state = 0
# state0 = 0


#print(state,state+shadow(2,H2_op,0))
# for i in tqdm(range(0,100000)):
#     #pool.apply_async(func=shadow(), args=(i, i + 1))
#     state = (state + shadow(6,H,4)/100000)
#     #state0 = state0 + shadow(6,H,0)/1000

# P_matrix = state.dot(state0)
# f1 = np.trace(P_matrix.dot(P_matrix))
# f2 = np.trace(P_matrix.dot(P_matrix.dot(P_matrix)))
# def f(x):
#     x=x[0]
#     return [x**3 - 3 * x * f1 + 2 * f2]

x0=[0.85]
#result = fsolve(f,x0)
#print('result:',result)
#H2_op_matrix = np.kron(Pauli(1),Pauli(1))
# energy1 = np.trace(state.dot(state.dot(state)))
# energy2 = np.trace(state.dot(state0))
# energy3 = np.trace(state0)
# energy4 = np.trace(H_matrix.dot(state0))
# energy5 = np.trace(H_matrix.dot(state))
# print(energy1,energy2,energy3,energy4,energy5,state0)
#print(energy1,energy2,energy3)
state = H_evolution(6,H,2)[-1]
print(state)
np.savetxt('t=2_state',state)
#0.0175762





# Create a Quantum Circuit acting on the q register


