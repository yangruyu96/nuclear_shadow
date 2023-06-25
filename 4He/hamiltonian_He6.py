import math

import numpy as np
import scipy as sc
from scipy import integrate
from sympy.physics.wigner import clebsch_gordan
from scipy.special import comb, perm
import itertools
from tqdm import trange
import time
from tqdm import tqdm
def Pauli(i):
    if i == 0:
        return np.array([[1,0],[0,1]])
    if i == 1:
        return np.array([[0,1],[1,0]])
    if i == 2:
        return np.array([[0,-1j],[1j,0]])
    if i == 3:
        return np.array([[1,0],[0,-1]])


def delta(a,b):
    if a==b:
        return 1
    else:
        return 0

def square(a,b):
    return 1/(math.sqrt(1+delta(a,b)))

def Jordan_creation(i):
    matrix_list = [] #1,2,0分别代表Z，产生算符，单位算符
    for k in range(0,i):
        matrix_list.append(1)
    matrix_list.append(2)
    for m in range(i+1,12):
        matrix_list.append(0)
    #print(matrix_list)
    # matrix = [[1]]
    # for k in range(0,i):
    #     matrix=np.kron(matrix,-Pauli(3))
    # matrix = np.kron(matrix,(Pauli(1) + 1j * Pauli(2))/2)
    # for m in range(i+1,12):
    #     matrix = np.kron(matrix,Pauli(0))
    return matrix_list

def Jordan_annilation(i):
    matrix_list = []  # 1,2,0分别代表Z，产生算符，单位算符,3代表消灭算符.
    for k in range(0, i):
        matrix_list.append(1)
    matrix_list.append(3)
    for m in range(i + 1, 12):
        matrix_list.append(0)
    #print(matrix_list)
    # matrix = [[1]]
    # for k in range(0,i):
    #     matrix=np.kron(matrix,-Pauli(3))
    # matrix = np.kron(matrix,(Pauli(1) - 1j * Pauli(2))/2)
    # for m in range(i+1,12):
    #     matrix = np.kron(matrix,Pauli(0))
    return matrix_list

#print(Jordan_annilation(12))

def Hamiltonian_Jordan(a,b,c,d):
    #print(a,b,c,d)
    matrix_a = Jordan_creation(a)
    matrix_b = Jordan_creation(b)
    matrix_c = Jordan_annilation(c)
    matrix_d = Jordan_annilation(d)
    return matrix_a + matrix_b + matrix_c + matrix_d
    #return np.dot(np.dot(np.dot(matrix_a,matrix_b),matrix_c),matrix_d)


def inner_V(x1,x2,x7,x8,J,M,T,M_T): # inner product of Hamiltonian interaction term
    x31,x41,x51,x61 = x1,x2,x7,x8 #case 1   3=1, 4=2, 5=7,  6=8
    x32, x42, x52, x62 = x2, x1, x7, x8  # case 2 3=2 4=1  5=7 6=8
    x33, x43, x53, x63 = x1, x2, x8, x7  # case 3 3=1, 4=2, 5=8,  6=7
    x34, x44, x54, x64 = x2, x1, x8, x7  # case 3
    j31,m31,u31 = x31[0],x31[1],x31[2]
    j32, m32, u32 = x32[0], x32[1], x32[2]
    j33, m33, u33 = x33[0], x33[1], x33[2]
    j34, m34, u34 = x34[0], x34[1], x34[2]
    j41, m41, u41 = x41[0], x41[1], x41[2]
    j42, m42, u42 = x42[0], x42[1], x42[2]
    j43, m43, u43 = x43[0], x43[1], x43[2]
    j44, m44, u44 = x44[0], x44[1], x44[2]
    j51, m51, u51 = x51[0], x51[1], x51[2]
    j52, m52, u52 = x52[0], x52[1], x52[2]
    j53, m53, u53 = x53[0], x53[1], x53[2]
    j54, m54, u54 = x54[0], x54[1], x54[2]
    j61, m61, u61 = x61[0], x61[1], x61[2]
    j62, m62, u62 = x62[0], x62[1], x62[2]
    j63, m63, u63 = x63[0], x63[1], x63[2]
    j64, m64, u64 = x64[0], x64[1], x64[2]
    if j31<= j41 and j61 <= j51:
        nor1 = 1
    else:
        nor1 = 0
    if j32<= j42 and j62 <= j52:
        nor2 = 1
    else:
        nor2 = 0
    if j33<= j43 and j63 <= j53:
        nor3 = 1
    else:
        nor3 = 0
    if j34<= j44 and j64 <= j54:
        nor4 = 1
    else:
        nor4 = 0
    # if nor1==1 and nor2 ==1 and nor3==1 and nor4==1:
    #     print(x1,x2,x3,x4)
    #print(nor1,nor2,nor3,nor4)
    coe1_left = clebsch_gordan(j31,j41,J,m31,m41,M) * clebsch_gordan(1/2,1/2,T,u31,u41,M_T)
    coe1_right = clebsch_gordan(j61,j51,J,m61,m51,M).conjugate() * clebsch_gordan(1/2,1/2,T,u61,u51,M_T).conjugate()
    #print(coe1_left,coe1_right,nor1,V_data(J,T,j31*2,j41*2,j61*2,j51*2))
    coe1 = coe1_left * coe1_right * nor1 * V_element(J,T,j31*2,j41*2,j61*2,j51*2)* square(j31,j41) * square(j51,j61)
    coe2_left = clebsch_gordan(j32, j42, J, m32, m42, M) * clebsch_gordan(1 / 2, 1 / 2, T, u32, u42, M_T)
    coe2_right = clebsch_gordan(j62, j52, J, m62, m52, M).conjugate() * clebsch_gordan(1 / 2, 1 / 2, T, u62, u52, M_T).conjugate()
    coe2 = coe2_left * coe2_right * nor2 * V_element(J,T,j32*2,j42*2,j62*2,j52*2)* (-1)* square(j32,j42) * square(j52,j62)
    coe3_left = clebsch_gordan(j33, j43, J, m33, m43, M) * clebsch_gordan(1 / 2, 1 / 2, T, u33, u43, M_T)
    coe3_right = clebsch_gordan(j63, j53, J, m63, m53, M).conjugate() * clebsch_gordan(1 / 2, 1 / 2, T, u63, u53, M_T).conjugate()
    coe3 = coe3_left * coe3_right * nor3 * V_element(J,T,j33*2,j43*2,j63*2,j53*2)* (-1)* square(j33,j43) * square(j53,j63)
    coe4_left = clebsch_gordan(j34, j44, J, m34, m44, M) * clebsch_gordan(1 / 2, 1 / 2, T, u34, u44, M_T)
    coe4_right = clebsch_gordan(j64, j54, J, m64, m54, M).conjugate() * clebsch_gordan(1 / 2, 1 / 2, T, u64, u54, M_T).conjugate()
    coe4 = coe4_left * coe4_right* nor4 * V_element(J,T,j34*2,j44*2,j64*2,j54*2)* square(j34,j44) * square(j54,j64)
    #print(V_element(J,T,j31*2,j41*2,j61*2,j51*2))
    #print(coe1,coe2,coe3,coe4,V_element(J,T,j32*2,j42*2,j62*2,j52*2))
    #print(V_element(J,T,j31*2,j41*2,j61*2,j51*2),V_element(J,T,j32*2,j42*2,j62*2,j52*2),V_element(J,T,j33*2,j43*2,j63*2,j53*2),V_element(J,T,j34*2,j44*2,j64*2,j54*2))
    # print(coe1_right,coe2_right,coe3_right,coe4_right)
    # print(j54,j64,m64,m54,J,M)
    # print(j34,j44,m34,m44,J,M)
    #print(clebsch_gordan(j54, j64, J, m64, m54,M))
    #print(clebsch_gordan(1 / 2, 1 / 2, T, u64, u54, M_T))
    return coe4 + coe3 + coe2 + coe1


def jordan_matrix(x):
    #print('x:',x)
    if x == 0:
        return Pauli(0)
    if x == 1:
        return Pauli(3)
    if x == 2:
        return (Pauli(1) + 1j * Pauli(2))/2
    if x == 3:
        return (Pauli(1) - 1j * Pauli(2))/2
    if x == 4:
        return -(Pauli(1) - 1j * Pauli(2))/2
    if x == 5:
        return -(Pauli(1) + 1j * Pauli(2))/2
    if x == 6:
        return np.array([[0,0],[0,0]])
    if x == 7:
        return np.array([[1,0],[0,0]])
    if x == 8:
        return np.array([[-1,0],[0,0]])
    if x==9:
        return np.array([[0,0],[0,1]])
    if x==10:
        return np.array([[0,0],[0,-1]])

def inner_V_sum(x1,x2,x7,x8,J,T):
    value = 0
    for M in range(-J,J+1):
        for M_T in range(-T,T+1):
            final = inner_V(x1, x2, x7, x8, J, M, T, M_T)
            value = value + final
            # if final !=0:
            #     print(J,M,T,M_T)
    return value

def inner_T(x1,x2,x7,x8,state):
    j,m,l = state[0],state[1],state[2]
    if state == x7:
        if x2 == x7 and x1 == x8:
            return 1
        if x1 == x7 and x2 == x8:
            return -1
        else:
            return 0
    if state == x8:
        if x1 == x8 and x2 == x7:
            return 1
        if x2 == x8 and x1 == x7:
            return -1
        else:
            return 0
    else:
        return 0

def V_data(J,T,a,b,c,d):
    #print(a,b,c,d)
    if J == 0 and T == 1:
        if a==3 and b==3 and c == 3 and d==3:
            return -2.74
        if a == 3 and b == 3 and c == 1 and d == 1:
            return -5.32
        if a == 1 and b == 1 and c == 1 and d == 1:
            return 0.34
        else:
            return 0
    if J == 1 and T == 1:
        if a == 3 and b == 1 and c == 3 and d == 1:
            return 0.86
        else:
            return 0
    if J == 2 and T == 1:
        if a == 3 and b == 3 and c == 3 and d == 3:
            return -0.65
        if a == 3 and b == 3 and c == 3 and d == 1:
            return -2.21
        if a == 3 and b == 1 and c == 3 and d == 1:
            return -1.14
        else:
            return 0
    if J == 1 and T == 0:
        if a == 3 and b == 3 and c == 3 and d == 3:
            return -3.14
        if a == 3 and b == 3 and c == 3 and d == 1:
            return 4.02
        if a == 3 and b == 3 and c == 1 and d == 1:
            return 1.09
        if a == 3 and b == 1 and c == 3 and d == 1:
            return -6.54
        if a == 3 and b == 1 and c == 1 and d == 1:
            return 1.39
        if a == 1 and b == 1 and c == 1 and d == 1:
            return -4.26
        else:
            return 0
    if J == 2 and T == 0:
        if a == 3 and b == 1 and c == 3 and d == 1:
            return -4.22
        else:
            return 0
    if J == 3 and T == 0:
        if a == 3 and b == 3 and c == 3 and d == 3:
            return -6.68
        else:
            return 0
    else:
        return 0

def j_angular(a):
    return a/2


def V_element(J,T,a,b,c,d):
    #return 0
    #print(V_data(J,T,b,a,d,c))
    if V_data(J,T,b,a,d,c) == 0:
        return V_data(J,T,d,c,b,a)#* ((-1)**((a+b+c+d)/2))
    return V_data(J,T,b,a,d,c)#* ((-1)**((a+b+c+d)/2))


def T_data(a):
    if a == 1:
        return 2.27
    if a == 3:
        return 1.63



def H_T_term(x1,x2,x7,x8):
    if x1 != x2 and x7 != x8:
        #print(x1,x2,x7,x8)
        list = [1,3]
        list_u = [1/2,-1/2]
        value = 0
        for a in list:
            ja = a/2
            for m in range(-a,a+1,2):
                ma = m/2
                for u in list_u:
                    state = [ja,ma,u]
                    value = value + T_data(a) * inner_T(x1,x2,x7,x8,state)
        return value
    else:
        return 0


def H_V_term(x1,x2,x7,x8):
    if x1 != x2 and x7 != x8:
        list = [1,3,5]
        list_u = [1/2,-1/2]
        value = 0
        for J in range(0,4):
            for T in range(0,2):
                x = inner_V_sum(x1, x2, x7, x8, J, T)
                value = value + x
                # if x != 0:
                #     print(x,J,T)
                #print(inner_V_sum(x1, x2, x7, x8, J, T))
        return value
    else:
        return 0

# x1 = [1/2,1/2,1/2]
# print(H_V_term(x1,x1,x1,x1))






#list = [1,3,5]
list = [1,3,5]
list_u = [-1/2]

def anti_index(x):
    if x==0:
        return [1/2,-1/2,-1/2]
    if x==1:
        return [1/2,-1/2,1/2]
    if x==2:
        return [1/2,1/2,-1/2]
    if x==3:
        return [1/2,1/2,1/2]
    if x==4:
        return [3/2,-3/2,-1/2]
    if x==5:
        return [3/2,-3/2,1/2]
    if x==6:
        return [3/2,-1/2,-1/2]
    if x==7:
        return [3/2,-1/2,1/2]
    if x==8:
        return [3/2,1/2,-1/2]
    if x==9:
        return [3/2,1/2,1/2]
    if x==10:
        return [3/2,3/2,-1/2]
    if x==11:
        return [3/2,3/2,1/2]




def jordan_dot(a,b):
    if a*b == 0:
        if a!= 0 :
            return a
        else:
            return b
    if a==1 and b==1:
        return 0
    if a==1 and b==2:
        return 2
    if a==1 and b==3:
        return 4
    if a==1 and b==4:
        return 3
    if a==1 and b==5:
        return 5
    if a==1 and b==6:
        return 6
    if a==1 and b==7:
        return 7
    if a==1 and b==8:
        return 8
    if a==1 and b==9:
        return 10
    if a==1 and b==10:
        return 9
    if a==2 and b==1:
        return 5
    if a==2 and b==2:
        return 6
    if a==2 and b==3:
        return 7
    if a==2 and b==4:
        return 8
    if a==2 and b==5:
        return 6
    if a==2 and b==6:
        return 6
    if a==2 and b==7:
        return 6
    if a==2 and b==8:
        return 6
    if a==2 and b==9:
        return 2
    if a==2 and b==10:
        return 5
    if a==3 and b==1:
        return 3
    if a==3 and b==2:
        return 9
    if a==3 and b==3:
        return 6
    if a==3 and b==4:
        return 6
    if a==3 and b==5:
        return 10
    if a==3 and b==6:
        return 6
    if a==3 and b==7:
        return 3
    if a==3 and b==8:
        return 4
    if a==3 and b==9:
        return 6
    if a==3 and b==10:
        return 6
    if a==4 and b==1:
        return 4
    if a==4 and b==2:
        return 10
    if a==4 and b==3:
        return 6
    if a==4 and b==4:
        return 6
    if a==4 and b==5:
        return 9
    if a==4 and b==6:
        return 6
    if a==4 and b==7:
        return 4
    if a==4 and b==3:
        return 3
    if a==4 and b==9:
        return 6
    if a==4 and b==10:
        return 6
    if a==5 and b==1:
        return 2
    if a==5 and b==2:
        return 6
    if a==5 and b==3:
        return 8
    if a==5 and b==4:
        return 7
    if a==5 and b==5:
        return 6
    if a==5 and b==6:
        return 6
    if a==5 and b==7:
        return 6
    if a==5 and b==8:
        return 6
    if a==5 and b==9:
        return 5
    if a==5 and b==10:
        return 2
    if a==6:
        return 6
    if a==7 and b==1:
        return 7
    if a==7 and b==2:
        return 2
    if a==7 and b==3:
        return 6
    if a==7 and b==4:
        return 6
    if a==7 and b==5:
        return 5
    if a==7 and b==6:
        return 6
    if a==7 and b==7:
        return 7
    if a==7 and b==8:
        return 8
    if a==7 and b==9:
        return 6
    if a==7 and b==10:
        return 6
    if a==8 and b==1:
        return 8
    if a==8 and b==2:
        return 5
    if a==8 and b==3:
        return 6
    if a==8 and b==4:
        return 6
    if a==8 and b==5:
        return 2
    if a==8 and b==6:
        return 6
    if a==8 and b==7:
        return 8
    if a==8 and b==8:
        return 7
    if a==8 and b==9:
        return 6
    if a==8 and b==10:
        return 6
    if a==9 and b==1:
        return 10
    if a==9 and b==2:
        return 6
    if a==9 and b==3:
        return 3
    if a==9 and b==4:
        return 4
    if a==9 and b==5:
        return 6
    if a==9 and b==6:
        return 6
    if a==9 and b==7:
        return 6
    if a==9 and b==8:
        return 6
    if a==9 and b==9:
        return 9
    if a==9 and b==10:
        return 10
    if a==10 and b==1:
        return 9
    if a==10 and b==2:
        return 6
    if a==10 and b==3:
        return 4
    if a==10 and b==4:
        return 3
    if a==10 and b==5:
        return 6
    if a==10 and b==6:
        return 6
    if a==10 and b==7:
        return 6
    if a==10 and b==8:
        return 6
    if a==10 and b==9:
        return 10
    if a==10 and b==10:
        return 9

# i_vector = anti_index(0)
# j_vector = anti_index(2)
# k_vector = anti_index(2)
# l_vector = anti_index(0)
#
# #print(H_V_term(i_vector,j_vector,k_vector,l_vector)/2)
# print(inner_V(i_vector,j_vector,k_vector,l_vector,0,0,1,-1))
H_T = []


for a in tqdm(range(0,12,1)):
    x = anti_index(a)
    #H =  (Jordan_creation(int(a/2)) + (Jordan_annilation(int(a/2))))
    H_T.append([T_data(x[0]*2),Jordan_creation(int(a)),(Jordan_annilation(int(a)))])
    #print(H_T)
#print(H_T[3])
# H_V = []
# print(inner_V(anti_index(0),anti_index(1),anti_index(3),anti_index(4),1,-1,0,0))
# print(inner_V(anti_index(4),anti_index(3),anti_index(1),anti_index(0),1,-1,0,0))
# print(H_V_term(anti_index(0),anti_index(1),anti_index(3),anti_index(4)))
# print(H_V_term(anti_index(4),anti_index(3),anti_index(1),anti_index(0)))


for a in tqdm(range(0,12,1)):
    for b in range(a+1,12,1):
        for d in range(0,12,1):
            for c in range(d+1,12,1):
                x1 = anti_index(a)
                x2 = anti_index(b)
                x4 = anti_index(d)
                x3 = anti_index(c)
                #print(len(H_T[0][-1]))
                H_T.append([complex(H_V_term(x1,x2,x3,x4)),Jordan_creation(int(a)),Jordan_creation(int(b)),Jordan_annilation(int(c)),Jordan_annilation(int(d))])

                #,Jordan_creation(int(a)),Jordan_creation(int(b)),Jordan_annilation(int(c)),Jordan_annilation(int(d))])
                #H = H + Hamiltonian_Jordan(int(a/2),int(b/2),int(c/2),int(d/2)) #+ H_T_term(x1,x2,x3,x4)
                #H.append(float(final ))9oiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiqw21111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111nqwassssssssssssssssssssssssssssssssssssssss
                #print(float(final))

state_left_list = []
state_right_list = []

for d in range(0,12,2):
    for c in range(d+2,12,2):
        for f in range(c+2,12,2):
                state_right_list.append([Jordan_creation(int(f)),Jordan_creation(int(c)),Jordan_creation(int(d))])

for a in tqdm(range(0,12,2)):
    for b in range(a+2,12,2):
        for e in range(b+2,12,2):
                state_left_list.append([Jordan_annilation(int(a)),Jordan_annilation(int(b)),Jordan_annilation(int(e))])

def matrix_element(list_left,list_h,list_right):
    left_len = len(list_left)
    h_len = len(list_h)
    right_len = len(list_right)
    matrix_list = []
    for v_every in range(0,len(list_h[-1])):
        matrix_every_pauli = 0
        for v_left in range(0,left_len):
            #print(list_left)
            #print(jordan_matrix(list_left[v_left][v_every]))
            #matrix_every_pauli = matrix_every_pauli.dot(jordan_matrix(list_left[v_left][v_every]))
            # if matrix_every_pauli == 0:
            #     continue
            matrix_every_pauli = jordan_dot(matrix_every_pauli,list_left[v_left][v_every])
        for v_h in range(1,h_len):
            #matrix_every_pauli = matrix_every_pauli.dot(jordan_matrix(list_h[v_h][v_every]))
        # if matrix_every_pauli == 0:
            #     continue
            matrix_every_pauli = jordan_dot(matrix_every_pauli, list_h[v_h][v_every])
        for v_right in range(0,right_len):
            #matrix_every_pauli = matrix_every_pauli.dot(jordan_matrix(list_right[v_right][v_every]))
            # if matrix_every_pauli == 0:
            #     continue
            matrix_every_pauli = jordan_dot(matrix_every_pauli, list_right[v_right][v_every])
        matrix_list.append(matrix_every_pauli)
    matrix_len = len(matrix_list)
    #initial_state = np.array([[0,1]])
    value = list_h[0]
    #if value != 0:
    for v_matrix in range(0,matrix_len):
        value = jordan_matrix(matrix_list[v_matrix])[1,1] * value
        if value == 0:
            #print(v_matrix)
            return 0
    return float(value)

len_left = len(state_left_list)
len_right = len(state_right_list)
H = []
#print(state_left_list)
for v_left in range(0,len_left):
    for v_right in tqdm(range(0,len_right)):
        H_term = 0
        for v_mid in range(0,len(H_T)):
            H_term = H_term + matrix_element(state_left_list[v_left],H_T[v_mid],state_right_list[v_right])
        H.append(H_term)

# #                 # if final != 0:
#                 #     print(a,b,c,d,final)
#                 #H.append(float(H_V_term(x1, x2, x3, x4) + H_T_term(x1, x2, x3, x4)))
# a,b,c,d=0,2,10,4
# x1 = anti_index(a)
# x2 = anti_index(b)
# x4 = anti_index(d)

# x3 = anti_index(c)
# final = inner_V(x1,x2,x3,x4,0,0,1,-1) - inner_V(x4,x3,x2,x1,0,0,1,-1).conjugate()
#
#
# print(final)
# for a in list:
#     for b in range(a,6,2):
#         ja = a/2
#         jb = b/2
#         #print(a,b)
#         for m1 in range(-a,a+1,2):
#             for m2 in range(m1,b+1,2):
#                 ma = m1/2
#                 mb = m2/2
#                 for ua in list_u:
#                     for ub in list_u:
#                         x1 = [ja,ma,ua]
#                         x2 = [jb,mb,ub]
#                         for d in list:
#                             for c in range(d,6,2):
#                                 jd = d / 2
#                                 jc = c / 2
#                                 for m4 in range(-d, d + 1, 2):
#                                     for m3 in range(m4, c + 1, 2):
#                                         md = m4 / 2
#                                         mc = m3 / 2
#                                         for ud in list_u:
#                                             for uc in list_u:
#                                                 x4 = [jd, md, ud]
#                                                 x3 = [jc, mc, uc]
#                                                 H.append(float(H_V_term(x1,x2,x3,x4) + H_T_term(x1,x2,x3,x4)))
#                                                 #print(x1,x2,x3,x4)

# print(H_T)
# H = np.array(H)
# l = int(len(H))
# l = math.sqrt(l)
# l = int(l)
# H = H.reshape(l,l)
# print("l:",l)
# print(type(H[1][1]))
# print(H)
# print(H_V)
np.savetxt('matrix_Hamiltonian_He7',H,fmt = '%s')
# print(H.shape)
#
# print(np.linalg.eig(H)[0])


#print(max(np.linalg.eig(H)[0]))

