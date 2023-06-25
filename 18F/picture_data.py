import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(p,x,y):
    k,b = p
    return k*x - y + b

p0 = [1/4,0]





fig = plt.figure('F18')
data = np.load('data_F18.npy')
x = data[0]
y = data[1]
print(1/y[-1])
x_linear = []
y_linear = []

para = leastsq(fun,p0,args=(x,y))
#print(para[0][0])
print(para)

for i in range(0,100):
    y_linear.append(i*para[0][0] + para[0][1])

plt.scatter(x,y,label = '${}^{18}$F_data')
plt.plot(y_linear,'r--',label = 'y=0.024x + 0.164')
plt.xlabel(r'shots $N/1000$')
plt.ylabel(r'accuracy $1/|\epsilon|$')
plt.legend()
plt.show()