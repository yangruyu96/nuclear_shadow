import math

import scipy
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from tqdm import tqdm
#
state_0 = np.zeros([2**10,2**10],dtype=complex)
state_0[0,0] = 1

state_11 = np.zeros([2**10,2**10],dtype=complex)
state_11[-1,-1] = 1


# state_1 = np.loadtxt('t=1_725',dtype=complex)
# for j in range(1,11):
#     state_x = np.loadtxt('t=1_{}_715'.format(j*10),dtype=complex)
#     state_1 = state_x + state_1
# state_1 = state_1/(np.trace(state_1))
#
# state_2 = np.loadtxt('t=2_725',dtype=complex)
# for j in range(1,11):
#     state_x = np.loadtxt('t=2_{}_715'.format(j*10),dtype=complex)
#     state_2 = state_x + state_2
# state_2 = state_2/(np.trace(state_2))


def state(i,n):
    state = np.loadtxt('shadow_{}_1k_Li6_{}_94_new'.format(n,i),dtype=complex) * 1000# + np.loadtxt('shadow_{}_100_Li6'.format(i),dtype=complex) * 100000 # + np.loadtxt('shadow_{}_200_Li6'.format(i),dtype=complex) * 100000\
   #+ np.loadtxt('shadow_{}_300_Li6'.format(i),dtype=complex) * 100000 + np.loadtxt('shadow_{}_400_Li6'.format(i),dtype=complex) * 100000
    #print(np.trace(state))
    return state/np.trace(state)

def state_ideal(i,n):
    state_vector = np.loadtxt('time_{}_Li6_new'.format(i),dtype=complex)
    state_vector_list = []
    #print(state_vector)
    for k in range(0, len(state_vector)):
        state_vector_list.append(state_vector[k])
    state_vector = np.array([state_vector_list])
    state_density_0 = state_vector.conjugate().T.dot(state_vector)
    return state_density_0

# def state(i):
#     state_vector = np.loadtxt('time_{}'.format(i),dtype=complex)
#     state_vector_list = []
#     #print(state_vector)
#     for k in range(0, len(state_vector)):
#         state_vector_list.append(state_vector[k])
#     state_vector = np.array([state_vector_list])
#     state_density_0 = state_vector.T.dot(state_vector.conjugate())
#     return state_density_0

def multi_trace(list):
    A = list[0]
    for i in range(1,len(list)):
        A = A.dot(list[i])
    return np.trace(A)

#H = np.loadtxt('hamiltonian_18F',dtype=complex)

h = np.loadtxt('matrix_Hamiltonian_Li6',dtype=complex)
L = len(h)
print('l:',np.sqrt(L))
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([64,64],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]

print(H)
np.savetxt('Li6_matrix',H)



# x_value = []
# y_value = []

# for scale in range(5,6):
#     for n in tqdm(range(0,100,10)):
#         state_list = []
#         x_value.append(((n+1)**(1) ))
#         for i in range(1,scale):
#             state_list.append(state(i,n))
#         for i in range(1, scale):
#             for j in range(1, scale):
#                 #if i!= j:
#                     state_list.append(state(i, n).dot(state(j, n)))
#         H_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
#         N_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
#         for i in range(0, len(state_list)):
#             for j in range(0, len(state_list)):
#                 H_direct[i, j] = multi_trace([state_list[i].T.conjugate(), H, state_list[j]])
#                 N_direct[i, j] = multi_trace([state_list[i].T.conjugate(), state_list[j]])
#         #print(H_direct.shape, N_direct.shape)
#         a, b = scipy.linalg.eig(H_direct, N_direct)
#         print(min(a))
#         a_ideal, b_ideal = scipy.linalg.eig(H)
#         y_ideal = min(a_ideal)
#         print(y_ideal)
#         #print(abs(min(a) - y_ideal))
#         y_value.append((abs(1/(min(a) - y_ideal))))


# d, v = scipy.linalg.eig(N_direct)
# print(d)
# a,b = scipy.linalg.eig(H_direct,N_direct)
# #print(a)
# d = np.diag(d)
# print(v.shape)
#v = np.delete(v, [-1,-2,-3,-4,-5,-6,-7,-8], axis=1)

# print(v.shape)
#


# print()

# N_matrix_th = v.T.conjugate().dot(N_direct.dot(v))
# H_matrix_th = v.T.conjugate().dot(H_direct.dot(v))
# a, b = scipy.linalg.eig(H_matrix_th, N_matrix_th)
# #print(a)
# print(min(a))
# a_ideal,b_ideal = scipy.linalg.eig(H)
# print(min(a_ideal))
# #y_value.append(abs(1/(min(a_ideal)-min(a))))
# # print(H_matrix)


# plt.scatter(x_value,y_value)
# plt.show()

# data = np.array([x_value,y_value])
# np.save('data_Li6',data)