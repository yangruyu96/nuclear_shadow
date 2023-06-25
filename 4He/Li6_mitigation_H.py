import scipy
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def state(i):
    state = np.loadtxt('shadow_0_1k_Li6_{}_94'.format(i),dtype=complex) #+ np.loadtxt('shadow_{}_200_Li6'.format(i),dtype=complex) + np.loadtxt('shadow_{}_300_Li6'.format(i),dtype=complex)\
   # +np.loadtxt('shadow_{}_400_Li6'.format(i),dtype=complex) \
     #      + np.loadtxt('shadow_{}_1_Li6'.format(i), dtype=complex)
    return state/np.trace(state)


def multi_trace(list):
    A = list[0]
    for i in range(1,len(list)):
        A = A.dot(list[i])
    return np.trace(A)

def state_ideal(i):
    state_vector = np.loadtxt('time_{}_Li6'.format(i),dtype=complex)
    state_vector_list = []
    #print(state_vector)
    for k in range(0, len(state_vector)):
        state_vector_list.append(state_vector[k])
    state_vector = np.array([state_vector_list])
    state_density_0 = state_vector.conjugate().T.dot(state_vector)
    return state_density_0

a = 2
b = 1

h = np.loadtxt('matrix_Hamiltonian_Li6',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([64,64],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]


print(multi_trace([state(a),state_ideal(b)]))
print(multi_trace([state_ideal(a),state(b)]))
print(multi_trace([state(a) - state_ideal(a),state(b) - state_ideal(b)]))


print(multi_trace([state(a),H,state(b)]))
print(multi_trace([state_ideal(a),H,state_ideal(b)]))
ga = multi_trace([state_ideal(a),state_ideal(b)])

#
#
# print(multi_trace([state(1),state(2),state(1),state(2),state(1),state(2)]))
# print(multi_trace([state_ideal(1),state_ideal(2),state_ideal(1),state_ideal(2),state_ideal(1),state_ideal(2)]))





x2 = multi_trace([state(a),state(b),state(a),state(b)])
x3 = multi_trace([state(a),state(b),state(a),state(b),state(a),state(b)])


def function(x):
    return (3/2) * x * x2 - x3 - (1/2)*x**3

#aa = fsolve(function,ga.real)
#print(aa)

#print(function(aa))

x_value_real = []
x_value_ima = []
y_value = []
x_value = []
for i in range(-1000,1000,1):
    for j in range(0,1000,1):
        x_value_real.append(i/1000)
        x_value_ima.append(j/1000)
        y_value.append(function((i/1000) + (j/1000)*1j))
        x_value.append(i/1000 + j/1000 * 1j)
print(abs(function(ga)))

plt.plot(x_value,y_value)
plt.show()

for i in range(0,len(y_value)):
    y_value[i] = abs(y_value[i])

print(min(y_value))

index = np.argmin(y_value)
print(index)
print(x_value[index])