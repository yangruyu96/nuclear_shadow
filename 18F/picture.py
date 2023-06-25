import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(p,x,y):
    k,b = p
    return k*x - y + b

p0 = [1/4,0]





fig = plt.figure('F18_diff')
data = np.load('data_F18.npy')
x = data[0]-1
y = data[1]

x_ss = x[3:12]
y_ss = y[3:12]

x_linear = []
y_linear = []

para = leastsq(fun,p0,args=(x_ss,y_ss))
print(para[0])
#print(para)
for i in range(4,11):
    x_linear.append(i)
    y_linear.append(i*para[0][0] + para[0][1])


plt.plot(x_linear,y_linear,'r--',label = 'y=0.384x-627')
plt.scatter(x,y,label = '${}^{18}$F_data')
plt.axvline(x=3.5,linestyle = '--',color = 'green')

plt.xlabel(r'number of evolved states',fontsize = 20)
plt.ylabel(r'accuracy $log(1/|\epsilon|)$',fontsize = 20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=18)
plt.show()