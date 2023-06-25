import scipy
import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
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


# def state(i):
#     state = np.loadtxt('shadow_{}'.format(i),dtype=complex)
#     return state/np.trace(state)

def state(i, n):
    state = np.loadtxt('shadow_{}_1k_F18_{}_94_new'.format(n, i),
                       dtype=complex) * 100  # + np.loadtxt('shadow_{}_1_He6'.format(i),dtype=complex) * 10000
    return state / np.trace(state)

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

h = np.loadtxt('hamiltonian_18F',dtype=complex)
print(h.shape)
H = np.zeros([256,256],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]

state_list = []

# for i in range(1,5):
#     state_list.append(state(i))

# for i in range(1,5):
#     for j in range(1,5):
#         state_list.append(state(i).dot(state(j)))

# for i in range(8,12):
#     for j in range(8,12):
#         for k in range(8,12):
#             state_list.append(state(i).dot(state(j).dot(state(k))))

H_direct = np.zeros([len(state_list),len(state_list)],dtype=complex)
N_direct = np.zeros([len(state_list),len(state_list)],dtype=complex)

# for i in range(0,len(state_list)):
#     for j in range(0,len(state_list)):
#         H_direct[i,j] = multi_trace([state_list[i].T.conjugate(),H,state_list[j]])
#         N_direct[i, j] = multi_trace([state_list[i].T.conjugate(), state_list[j]])

x_value = []
y_value = []
for scale in range(5,6):
    for n in range(0,100,10):
        state_list = []
        x_value.append((n+1)**(1))
        for i in range(1,scale):
            state_list.append(state(i,n))
        for i in range(1, scale):
            for j in range(1, scale):
                state_list.append(state(i, n).dot(state(j, n)))
            #print(state_list)
        H_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
        N_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
        for i in range(0, len(state_list)):
            for j in range(0, len(state_list)):
                H_direct[i, j] = multi_trace([state_list[i].T.conjugate(), H, state_list[j]])
                N_direct[i, j] = multi_trace([state_list[i].T.conjugate(), state_list[j]])
        #print(H_direct.shape, N_direct.shape)
        a, b = scipy.linalg.eig(H_direct, N_direct)
        #print(min(a))
        a_ideal, b_ideal = scipy.linalg.eig(H)
        y_ideal = min(a_ideal)
        #print(min(a))
        #print(y_ideal)
        #print(abs(min(a) - y_ideal))
        y_value.append(abs(1/(min(a) - y_ideal)) )


plt.scatter(x_value,y_value)

plt.show()
data = np.array([x_value,y_value])
#np.save('data_F18',data)

#print(H_direct.shape,N_direct.shape)

# a,b = scipy.linalg.eig(H_direct,N_direct)
# #np.savetxt('H_direct_4-1',H_direct)
# #np.savetxt('N_direct_4-1',N_direct)
#
# print(min(a))
# a,b = scipy.linalg.eig(H)
# print(min(a))


