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
    state =  np.loadtxt('shadow_{}_1k_Li8_{}_94_new'.format(n,i),dtype=complex) #+ np.loadtxt('shadow_0_1w_Li8_{}'.format(i),dtype=complex)
    #+ np.loadtxt('shadow_{}_30w_Li8'.format(i),dtype=complex) + np.loadtxt('shadow_{}_40w_Li8'.format(i),dtype=complex)\
    #+ np.loadtxt('shadow_{}_60w_Li8'.format(i),dtype=complex)  + np.loadtxt('shadow_{}_70w_Li8'.format(i),dtype=complex)
    return state/np.trace(state)

def state_ideal(i,n):
    state_vector = np.loadtxt('time_{}_Li8'.format(i),dtype=complex)
    state_vector_list = []
    #print(state_vector)
    for k in range(0, len(state_vector)):
        state_vector_list.append(state_vector[k])
    state_vector = np.array([state_vector_list])
    state_density_0 = state_vector.conjugate().T.dot(state_vector)
    return state_density_0


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



refe = 2
def overlap_basis(n,i,j=refe):
    state_1 = state(i,n)
    state_2 = state(j,n)
    x = multi_trace([state_1,state_2])
    #print(x)
    # state_1 = state_vector(i)
    # state_2 = state_vector(j)
    # x = state_1.conjugate().dot(state_2.T)
    return math.sqrt(x)

def overlap_basis_ideal(n,i,j=refe):
    state_1 = state_ideal(i,n)
    state_2 = state_ideal(j,n)
    x = multi_trace([state_1,state_2])
    #print(x)
    # state_1 = state_vector(i)
    # state_2 = state_vector(j)
    # x = state_1.conjugate().dot(state_2.T)
    return math.sqrt(x)

# def overlap(i,j):
#     return state_vector(i).dot(state_vector(j).T.conjugate())

def overlap(i,j,n):
    if i != refe and j != refe:
        state_1 = state(i,n)
        state_2 = state(j,n)
        state_10 = state(refe,n)
        x = multi_trace([state_1,state_2,state_10])
        y1 = overlap_basis(n,i,refe).conjugate()
        y2 = overlap_basis(n,j,refe)
        return (x/(y1 * y2))
    else:
        if i != refe:
            return overlap_basis(n,i,refe)
        else:
            if j != refe:
                return overlap_basis(n,j,refe).conjugate()
            else:
                return (multi_trace([state(refe,n),state(7,n),state(refe,n)])/multi_trace([state(refe,n),state(7,n)]))

def overlap_ideal(i,j,n):
    if i != refe and j != refe:
        state_1 = state_ideal(i,n)
        state_2 = state_ideal(j,n)
        state_10 = state_ideal(refe,n)
        x = multi_trace([state_1,state_2,state_10])
        y1 = overlap_basis_ideal(n,i,refe).conjugate()
        y2 = overlap_basis_ideal(n,j,refe)
        return (x/(y1 * y2))
    else:
        if i != refe:
            return overlap_basis_ideal(n,i,refe)
        else:
            if j != refe:
                return overlap_basis_ideal(n,j,refe).conjugate()
            else:
                return (multi_trace([state_ideal(refe,n),state_ideal(7,n),state_ideal(refe,n)])/multi_trace([state_ideal(refe,n),state_ideal(7,n)]))



h = np.loadtxt('matrix_Hamiltonian_Li8',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([128,128],dtype=complex)
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
        if i != refe:
            x = multi_trace([H,state(i,n),state(refe,n),state(j,n)])
            y = multi_trace([state(i,n),state(refe,n)])
            return x/y
        else:
            x = multi_trace([H,state(i,n),state(7,n),state(j,n)])
            y = multi_trace([state(i,n),state(7,n)])
            return x/y

def H_overlap_ideal(j,i,n):
    if i != j:
        x = multi_trace([H,state_ideal(i,n),state_ideal(j,n)])
        y = overlap_ideal(i,j,n)
        return x/y
    else:
        if i != refe:
            x = multi_trace([H,state_ideal(i,n),state_ideal(refe,n),state_ideal(j,n)])
            y = multi_trace([state_ideal(i,n),state_ideal(refe,n)])
            return x/y
        else:
            x = multi_trace([H,state_ideal(i,n),state_ideal(7,n),state_ideal(j,n)])
            y = multi_trace([state_ideal(i,n),state_ideal(7,n)])
            return x/y


#print(H_overlap(1,1) * overlap(1,1))

H_matrix = np.zeros([7,7],dtype=complex)
N_matrix = np.zeros([7,7],dtype=complex)

H_matrix_ideal  = np.zeros([12,12],dtype=complex)
N_matrix_ideal = np.zeros([12,12],dtype=complex)

# print(overlap(3,4))
# x = state_vector(3).T.conjugate().dot(state_vector(4))
# print(x)
# print(multi_trace([state(3),state(4)]))


#print(overlap(10,10) * overlap(10,10))
x_value = []
y_value = []
for n in tqdm(range(360,370,10)):
    x_value.append(((n + 1) ** (1)))
    for i in tqdm(range(0, 7)):
        for j in tqdm(range(0, 7)):
            H_matrix[i, j] = H_overlap(i + 1, j + 1,n)
            N_matrix[i, j] = overlap(i + 1, j + 1,n)
            H_matrix_ideal[i, j] = H_overlap_ideal(i + 1, j + 1, n)
            N_matrix_ideal[i, j] = overlap_ideal(i + 1, j + 1, n)
    #a, b = scipy.linalg.eig(H_matrix, N_matrix)
    d, v = scipy.linalg.eig(H_matrix)
    print(d)
    a,b = scipy.linalg.eig(H_matrix,N_matrix)
    print(a)
    d = np.diag(d)
    # print(v.shape)
    #v = np.delete(v, [-1,-2,-3,-4,-5,-6], axis=1)

    # print(v.shape)
    #
    N = np.delete(N_matrix, [-1], axis=0)
    N = np.delete(N, [-1], axis=1)

    # print()

    N_matrix_th = v.T.conjugate().dot(N_matrix.dot(v))
    H_matrix_th = v.T.conjugate().dot(H_matrix.dot(v))
    a, b = scipy.linalg.eig(H_matrix_th, N_matrix_th)
    #print(a)
    print(min(a))
    a_ideal,b_ideal = scipy.linalg.eig(H)
    y_value.append(abs(1/(min(a_ideal)-min(a))))
# print(H_matrix)
# print(N_matrix)
#
#np.savetxt('H_18O_matrix_10_87-1',H_matrix)
# #np.savetxt('N_18O_matrix_10_87-1',N_matrix)
plt.scatter(x_value,y_value)
plt.show()

# a,b = scipy.linalg.eig(H_matrix,N_matrix)
# print(min(a))
#
# a,b = scipy.linalg.eig(H_matrix_ideal,N_matrix_ideal)
# print(min(a))
#
# print(H_matrix_ideal - H_matrix)
# #a,b = scipy.linalg.eig(H)
# a,b = scipy.linalg.eig(N_matrix)
# print(a)
#print(a)

# d,v = scipy.linalg.eig(N_matrix)
# print(d)
# d = np.diag(d)
# #print(v.shape)
# v = np.delete(v,[-1,-2,-3,-4,-5,-6],axis=1)
#
#
# #print(v.shape)
# #
# N = np.delete(N_matrix,[-1],axis=0)
# N = np.delete(N,[-1],axis=1)
#
# #print()
#
# N_matrix_th = v.T.conjugate().dot(N_matrix.dot(v))
#
# H_matrix_th = v.T.conjugate().dot(H_matrix.dot(v))
#
# a,b = scipy.linalg.eig(H_matrix_th,N_matrix_th)
# print(a)
#print(N_matrix_th)