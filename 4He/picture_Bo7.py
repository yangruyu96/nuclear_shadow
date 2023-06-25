import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(k,x,y):
    return k*x - y

p0 = [1/4]





fig = plt.figure('Bo7')
data = np.load('data_Bo7.npy')
x = data[0]
y = data[1]
print(1/y[-1])
x_linear = []
y_linear = []

para = leastsq(fun,p0,args=(x,y))
print(para[0][0])

for i in range(0,100):
    y_linear.append(i*para[0][0])


plt.plot(y_linear,'r--',label = 'y=0.148x')
plt.scatter(x,y,label = '${}^7$Bo_data')


plt.xlabel('shots $N/100$',fontsize=20)
plt.ylabel('accuracy $1/|\epsilon|$',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=18)

plt.show()