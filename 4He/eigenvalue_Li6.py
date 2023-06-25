import math
import scipy
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm
import matplotlib.pyplot as plt
# def state(i):
#     state_vector = np.loadtxt('time_{}'.format(i),dtype=complex)
#     state_vector_list = []
#     #print(state_vector)
#     for k in range(0, len(state_vector)):
#         state_vector_list.append(state_vector[k])
#     state_vector = np.array([state_vector_list])
#     state_density_0 = state_vector.T.dot(state_vector.conjugate())
#     return state_density_


def state(i,n):
    state = np.loadtxt('shadow_{}_1k_Li6_{}_94_new'.format(n,i),dtype=complex)*1000 #+ np.loadtxt('shadow_{}_200_Li6'.format(i),dtype=complex) * 100000 + np.loadtxt('shadow_{}_300_Li6'.format(i),dtype=complex) * 100000\
    #+np.loadtxt('shadow_{}_400_Li6'.format(i),dtype=complex) \
     #      + np.loadtxt('shadow_{}_1_Li6'.format(i), dtype=complex)
    return state

def state_vector(i):
    x = np.loadtxt('time_{}_Li6'.format(i),dtype=complex)
    return x

# def state(i):
#     state = np.loadtxt('t={}_725'.format(i),dtype=complex)
#     return state/(np.trace(state))

def multi_trace(list):
    A = list[0]
    for i in range(1,len(list)):
        A = A.dot(list[i])
    return np.trace(A)



def overlap_basis(n,i,j=4):
    state_1 = state(i,n)
    state_2 = state(j,n)
    x = multi_trace([state_1,state_2])
    # state_1 = state_vector(i)
    # state_2 = state_vector(j)
    # x = state_1.conjugate().dot(state_2.T)
    return math.sqrt(x)

# def overlap(i,j):
#     return state_vector(i).dot(state_vector(j).T.conjugate())

def overlap(i,j,n):
    if i != 4 and j != 4:
        state_1 = state(i,n)
        state_2 = state(j,n)
        state_10 = state(4,n)
        x = multi_trace([state_1,state_2,state_10])
        y1 = overlap_basis(n,i,4).conjugate()
        y2 = overlap_basis(n,j,4)
        return (x/(y1 * y2))
    else:
        if i != 4:
            return overlap_basis(n,i,4)
        else:
            if j != 4:
                return overlap_basis(n,j,4).conjugate()
            else:
                return (multi_trace([state(4,n),state(3,n),state(4,n)])/multi_trace([state(4,n),state(3,n)]))


h = np.loadtxt('matrix_Hamiltonian_Li6',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([64,64],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]

#print(overlap(5,2)*overlap(2,5))

def H_overlap(j,i,n):
    if i != j:
        x = multi_trace([H,state(i,n),state(j,n)])
        y = overlap(i,j,n)
        return x/y
    else:
        if i != 4:
            x = multi_trace([H,state(i,n),state(4,n),state(j,n)])
            y = multi_trace([state(i,n),state(4,n)])
            return x/y
        else:
            x = multi_trace([H,state(i,n),state(3,n),state(j,n)])
            y = multi_trace([state(i,n),state(3,n)])
            return x/y

#print(H_overlap(1,1) * overlap(1,1))

H_matrix = np.zeros([4,4],dtype=complex)
N_matrix = np.zeros([4,4],dtype=complex)

# print(overlap(3,4))
# x = state_vector(3).T.conjugate().dot(state_vector(4))
# print(x)
# print(multi_trace([state(3),state(4)]))


#print(overlap(10,10) * overlap(10,10))

x_value = []
y_value = []
for n in tqdm(range(0,50,10)):
    x_value.append(((n + 1) ** (1)))
    for i in tqdm(range(0, 4)):
        for j in tqdm(range(0, 4)):
            H_matrix[i, j] = H_overlap(i + 1, j + 1,n)
            N_matrix[i, j] = overlap(i + 1, j + 1,n)
    a, b = scipy.linalg.eig(H_matrix, N_matrix)
    print(min(a))
    a_ideal,b_ideal = scipy.linalg.eig(H)
    y_value.append(abs(min(a_ideal)-min(a)))

# print(H_matrix)
# print(N_matrix)
#
#np.savetxt('H_18O_matrix_10_87-1',H_matrix)
#np.savetxt('N_18O_matrix_10_87-1',N_matrix)

plt.scatter(x_value,y_value)
plt.show()
