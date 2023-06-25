import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def fun(p,x,y):
    k,b =p
    return k*x - y + b

p0 = [1/4,0]





fig = plt.figure('Bo7_bar')
data = np.load('data_Bo7_bar.npy')
x = data[0]
y = data[1]
x_linear = []
y_linear = []


print(len(x))

# plt.scatter(x,y)
# plt.plot(y_linear,'r--')
# plt.show()



x_average = []
y_average = []
y_error = []
y_down = []
for i in range(0,10):
    x_reduce = []
    y_reduce = []
    for j in range(0,6):
        x_reduce.append(x[10*j + i])
        y_reduce.append(y[10*j + i])
    x_average.append(x[i])
    y_average.append(np.mean(y_reduce))
    y_error.append(np.std(y_reduce))
    y_down.append(y_average[-1] - y_error[-1])

x_average = np.array(x_average)
print(type(x_average))
print(type(x))

para = leastsq(fun,p0,args=(x,y))
print(para[0])

for i in range(0,100):
    y_linear.append(i*para[0][0] + para[0][1])
print(y_down)
#print(x_average)
p0 = [1/4,0]
#para = leastsq(fun,p0,args=(x,y))
para = leastsq(fun,p0,args=(x_average,y_down))
print(para[0])

y_linear_error=[]

for i in range(0,100):
    y_linear_error.append(i*para[0][0] + para[0][1])


plt.errorbar(x_average,y_average,yerr=y_error,fmt='ob',ecolor='k',capsize=3,label = '${}^7$Bo_data')
plt.plot(y_linear_error,'g--',label = '0.121x+0.323')
plt.plot(y_linear,'r--',label = 'y=0.154x-0.075')
plt.legend()
plt.xlabel('shots $N/100')
plt.ylabel('accuracy $|1/\epsilon$')
plt.show()

