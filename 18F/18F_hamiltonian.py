import math

import numpy as np

hamiltonian = np.loadtxt('matrix_Hamiltonian_18F_625',dtype=complex)

hamiltonian = hamiltonian.reshape([144,144])

a,b = np.linalg.eig(hamiltonian)

print(min(a))

np.savetxt('hamiltonian_18F',hamiltonian)

#print(math.sqrt(20736))