import numpy as np

state1 = 0
for i in range(10,101,10):
    try:
        state = np.loadtxt('t=2_{}'.format(i),dtype=complex)
        state1 = state1 + state
    except:
        state1 = state1
state1 = state1/10
state2 = 0
for i in range(10,101,10):
    try:
        state = np.loadtxt('t=4_{}'.format(i),dtype=complex)
        state2 = state2 + state
    except:
        state2 = state2
state2 = state2/10
state1_vector = np.loadtxt('t=2_state',dtype=complex)

state1_vector_list = []
for i in range(0,len(state1_vector)):
    state1_vector_list.append(state1_vector[i])
state1_vector = np.array([state1_vector_list])
state4_vector = np.loadtxt('t=4_state',dtype=complex)
state4_vector_list = []
for i in range(0,len(state4_vector)):
    state4_vector_list.append(state4_vector[i])
state4_vector = np.array([state4_vector_list])
#print(state4_vector.shape)
#print(state1_vector.shape)
#print(state1_vector.dot(state1_vector.T.conjugate()))
#print(state1_vector.T.conjugate().dot(state1_vector))
state1_density = state1_vector.T.conjugate().dot(state1_vector)
state4_density = state4_vector.T.conjugate().dot(state4_vector)
#print(state1_density.shape)
# state1_vector = state1_vector.T.conjugate().dot(state1_vector)
H = np.loadtxt('test',dtype=complex)
# #print(np.trace(state1_vector.dot(H)))
# print(np.trace(state1.dot(H)))
# print(np.trace(state1_density.dot(H)))
print(np.trace(state1.dot(state4_density)))
print(np.trace(state1_density.dot(state4_density)))
print(np.trace(state1.dot(state2)))
# print(np.trace(state1.dot(state1_vector)))
# print(np.trace(state1_vector.dot(state1_vector)))
# print(state1_vector.shape)
print(np.trace(state1.dot(H.dot(state2))))
print(np.trace(state1_density.dot(H.dot(state4_density))))
#print(np.trace(state1))

