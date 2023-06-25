import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(p,x,y):
    k,b =p
    return k*x - y + b

p0 = [1/4,0]





fig = plt.figure('Li8_diff_sub')
data = np.load('data_Li8_diff_sub_523.npy')
x = data[0]-1
print(x)
#print(y)
x_ss = x[6:18]
print(x_ss)

#print(x[0])
#x = np.delete(x[0],obj=None)
#print(x)
y = data[1]
y_ss = y[6:18]
print(y_ss)
# for i in range(0,len(x)-1):
#     x[i] = x[i+1]
#     y[i] = y[i+1]
#y = np.delete(y[0],0)
x_linear = []
y_linear = []

para = leastsq(fun,p0,args=(x_ss,y_ss))
print(para)

for i in range(7,20):
    x_linear.append(i)
    y_linear.append(i*para[0][0] + para[0][1])


plt.plot(x_linear,y_linear,'r--',label = 'y=0.170x-0.640')
plt.scatter(x,y,label = '${}^8$Li_data')
plt.xticks([0,2,4,6,8,10,12,14,16,18])
plt.axvline(x=6.5,linestyle = '--',color = 'green')
#plt.xlim(7,14)
plt.xlabel('number of evolved states',fontsize = 20)
plt.ylabel('accuracy $log(1/|\epsilon|)$',fontsize = 20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.ylim(-2.3,2.5)
plt.legend(fontsize=18)
plt.show()