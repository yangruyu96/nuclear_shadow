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

print(multi_trace([state(a),state_ideal(b)]))
print(multi_trace([state_ideal(a),state(b)]))
print(multi_trace([state(a) - state_ideal(a),state(b) - state_ideal(b)]))


print(multi_trace([state(a),state(b)]))
print(multi_trace([state_ideal(a),state_ideal(b)]))
ga = multi_trace([state_ideal(a),state_ideal(b)])

#
#
# print(multi_trace([state(1),state(2),state(1),state(2),state(1),state(2)]))
# print(multi_trace([state_ideal(1),state_ideal(2),state_ideal(1),state_ideal(2),state_ideal(1),state_ideal(2)]))





x2 = multi_trace([state(a),state(b),state(a),state(b)]).real
x3 = multi_trace([state(a),state(b),state(a),state(b),state(a),state(b)]).real


def function(x):
    return (3/2) * x * x2 - x3 - (1/2)*x**3

aa = fsolve(function,ga.real)
print(aa)

print(function(aa))

x_value = []
y_value = []
for i in range(0,100,1):
    x_value.append(i/100)
    y_value.append(function(i/100))

plt.plot(x_value,y_value)
plt.show()