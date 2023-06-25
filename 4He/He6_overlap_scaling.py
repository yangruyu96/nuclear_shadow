import scipy
import numpy as np
from scipy.optimize import fsolve
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



def state(i,n):
    state = np.loadtxt('shadow_{}_1k_Li8_{}_94_new'.format(n,i),dtype=complex)
    return state/np.trace(state)
def state_ideal(i):
    state_vector = np.loadtxt('time_{}_He6'.format(i),dtype=complex)
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

h = np.loadtxt('matrix_Hamiltonian_He6',dtype=complex)
L = len(h)
h = np.array(h)
h = h.reshape([int(np.sqrt(L)),int(np.sqrt(L))])
H = np.zeros([16,16],dtype=complex)
for i in range(0,len(h)):
    for j in range(0,len(h)):
        H[i,j] = h[i,j]



# for i in range(1,4):
#     state_list.append(state(i))

x_value = []
y_value = []
for n in range(0,10,1):
    state_overlap = multi_trace([state(1,n),state(2,n)])
    state_overlap_ideal = multi_trace([state_ideal(1),state_ideal(2)])
    x_value.append((n + 1) ** (1))
    print(state_overlap)
    print(state_overlap_ideal)
    y_ideal = state_overlap_ideal
    y_value.append(abs(1/(state_overlap - y_ideal)))

plt.scatter(x_value,y_value)
plt.ylim(0,100)
plt.show()

# for i in range(8,12):
#     for j in range(8,12):
#         for k in range(8,12):
#             state_list.append(state(i).dot(state(j).dot(state(k))))


#np.savetxt('H_direct_4-1',H_direct)
#np.savetxt('N_direct_4-1',N_direct)


a,b = scipy.linalg.eig(H)
print(min(a))

