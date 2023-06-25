import math
import scipy
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm

# def state(i):
#     state_vector = np.loadtxt('time_{}'.format(i),dtype=complex)
#     state_vector_list = []
#     #print(state_vector)
#     for k in range(0, len(state_vector)):
#         state_vector_list.append(state_vector[k])
#     state_vector = np.array([state_vector_list])
#     state_density_0 = state_vector.T.dot(state_vector.conjugate())
#     return state_density_0


state_0 = np.zeros([2**10,2**10],dtype=complex)
state_0[0,0] = 1

state_11 = np.zeros([2**10,2**10],dtype=complex)
state_11[-1,-1] = 1
state_1 = 0
for j in range(1,2):
    state_x = np.loadtxt('t=1_{}'.format(j*10),dtype=complex)
    state_1 = state_1 + state_x

state_1 = state_1/(np.trace(state_1))

state_2 = 0
for j in range(1,2):
    state_x = np.loadtxt('t=2_{}'.format(j*10),dtype=complex)
    state_2 = state_2 + state_x

state_2 = state_2/(np.trace(state_2))


state_3 = 0
for j in range(1,2):
    state_x = np.loadtxt('t=3_{}'.format(j*10),dtype=complex)
    state_3 = state_3 + state_x

state_3 = state_3/(np.trace(state_3))

state_4 = 0
for j in range(1,2):
    state_x = np.loadtxt('t=4_{}'.format(j*10),dtype=complex)
    state_4 = state_4 + state_x

state_4 = state_4/(np.trace(state_4))





def state(i):
    if i == 1:
        return state_1
    if i == 2:
        return state_2
    if i == 3:
        return state_3
    if i == 4:
        return state_4

def state_vector(i):
    x = np.loadtxt('time_{}'.format(i),dtype=complex)
    return x

# def state(i):
#     state = np.loadtxt('t={}_725'.format(i),dtype=complex)
#     return state/(np.trace(state))

def multi_trace(list):
    A = list[0]
    for i in range(1,len(list)):
        A = A.dot(list[i])
    return np.trace(A)



def overlap_basis(i,j=4):
    state_1 = state(i)
    state_2 = state(j)
    x = multi_trace([state_1,state_2])
    # state_1 = state_vector(i)
    # state_2 = state_vector(j)
    # x = state_1.conjugate().dot(state_2.T)
    return math.sqrt(x)

# def overlap(i,j):
#     return state_vector(i).dot(state_vector(j).T.conjugate())

def overlap(i,j):
    if i != 4 and j != 4:
        state_1 = state(i)
        state_2 = state(j)
        state_10 = state(4)
        x = multi_trace([state_1,state_2,state_10])
        y1 = overlap_basis(i,4).conjugate()
        y2 = overlap_basis(j,4)
        return (x/(y1 * y2))
    else:
        if i != 4:
            return overlap_basis(i,4)
        else:
            if j != 4:
                return overlap_basis(j,4).conjugate()
            else:
                return (multi_trace([state(4),state(3),state(4)])/multi_trace([state(4),state(3)]))


H = np.loadtxt('test',dtype=complex)

#print(overlap(5,2)*overlap(2,5))

def H_overlap(j,i):
    if i != j:
        x = multi_trace([H,state(i),state(j)])
        y = overlap(i,j)
        return x/y
    else:
        if i != 4:
            x = multi_trace([H,state(i),state(4),state(j)])
            y = multi_trace([state(i),state(4)])
            return x/y
        else:
            x = multi_trace([H,state(i),state(3),state(j)])
            y = multi_trace([state(i),state(3)])
            return x/y

#print(H_overlap(1,1) * overlap(1,1))

H_matrix = np.zeros([4,4],dtype=complex)
N_matrix = np.zeros([4,4],dtype=complex)

# print(overlap(3,4))
# x = state_vector(3).T.conjugate().dot(state_vector(4))
# print(x)
# print(multi_trace([state(3),state(4)]))


#print(overlap(10,10) * overlap(10,10))
for i in tqdm(range(0,4)):
    for j in tqdm(range(0,4)):
        H_matrix[i,j] = H_overlap(i+1,j+1)
        N_matrix[i,j] = overlap(i+1,j+1)

# print(H_matrix)
# print(N_matrix)
#
np.savetxt('H_18O_matrix_10_87-1',H_matrix)
np.savetxt('N_18O_matrix_10_87-1',N_matrix)


a,b = scipy.linalg.eig(H_matrix,N_matrix)
print(min(a))