import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(p,x,y):
    k,b =p
    return k*x - y + b

p0 = [1/4,0]





fig = plt.figure('Li8')
data = np.load('data_Li8_523.npy')
x = data[0]
y = data[1]
print(1/y[-1])
x_linear = []
y_linear = []

para = leastsq(fun,p0,args=(x,y))
print(para)

for i in range(0,300):
    y_linear.append(i*para[0][0] + para[0][1])


plt.plot(y_linear,'r--',label='y=0.0044x+0.108')
plt.scatter(x,y,label = '${}^8$Li_data')


plt.xlabel('shots $N/1000$',fontsize=20)
plt.ylabel('accuracy $1/|\epsilon|$',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=18)
plt.show()