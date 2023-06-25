import math

import scipy
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#
state_0 = np.zeros([2**5,2**5],dtype=complex)
state_0[0,0] = 1

state_11 = np.zeros([2**5,2**5],dtype=complex)
state_11[1,1] = 1


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



def state(i,n,bar):
    if i == 0:
        return state_0
    else:
        state = np.loadtxt('shadow_{}_1h_Bo7_{}_94_new_bar_{}'.format(n,i,bar),dtype=complex) * 100 # + np.loadtxt('shadow_{}_1_He6'.format(i),dtype=complex) * 10000
        return state/np.trace(state)



# def state(i,n):
#     state_vector = np.loadtxt('time_{}_Bo7'.format(i),dtype=complex)
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

h = np.loadtxt('matrix_Hamiltonian_Bo7',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([32,32],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]



# for i in range(1,4):
#     state_list.append(state(i))

x_value = []
y_value = []
#scale = 5
for scale in range(90,100):
    for bar in range(1,2):
        for n in range(10,20,10):
            state_list = []
            x_value.append(scale)
            # for i in range(1,scale):
            #     state_list.append(state(i,n,bar))
            for i in range(1, scale):
                for j in range(1, scale):
                    state_list.append(state(i, n,bar).dot(state(j, n,bar)))
                #print(state_list)
            H_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
            N_direct = np.zeros([len(state_list), len(state_list)], dtype=complex)
            for i in range(0, len(state_list)):
                for j in range(0, len(state_list)):
                    H_direct[i, j] = multi_trace([state_list[i].T.conjugate(), H, state_list[j]])
                    N_direct[i, j] = multi_trace([state_list[i].T.conjugate(), state_list[j]])
            #print(H_direct.shape, N_direct.shape)
            N_pinv = np.linalg.pinv(N_direct)
            a,b = scipy.linalg.eig(N_pinv @ H_direct)
        #a, b = scipy.linalg.eig(H_direct, N_direct)
            print(scale,min(a))
            #print(min(a))
            a_ideal, b_ideal = scipy.linalg.eig(H)
            y_ideal = min(a_ideal)
            print(min(a))
            print(y_ideal)
            #print(abs(min(a) - y_ideal))
            y_value.append(math.log10(abs(1/(min(a) - y_ideal) )))
            del state_list

plt.scatter(x_value,y_value)
plt.show()

# for i in range(8,12):
#     for j in range(8,12):
#         for k in range(8,12):
#             state_list.append(state(i).dot(state(j).dot(state(k))))


#np.savetxt('H_direct_4-1',H_direct)
#np.savetxt('N_direct_4-1',N_direct)


a,b = scipy.linalg.eig(H)
print(min(a))
#print(multi_trace([H,state_0]))

data = np.array([x_value,y_value])
np.save('data_Bo7_scale_90——100',data)

#print(multi_trace([H,state(1,7000)]))
#print(multi_trace([H,state(3,7000)]))