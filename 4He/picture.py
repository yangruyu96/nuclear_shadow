# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
#
# #fig= plt.figure('abc')
#
# data = np.load('data_sub_diag_Bo7.npy')
# print(data[1])
# draw_circle = plt.Circle((0.5, 0.5), 0.3,fill=False)
#
# x = data[0]
# y = data[1]
# plt.scatter(x,y,marker = '.',s=10)
# #plt.ylim(-1,0)
#
# plt.gcf().gca().add_artist(draw_circle)
# #axes.set_aspect(1)
# plt.show()
#
#
# figure, axes = plt.subplots()
# draw_circle = plt.Circle((0.5, 0.5), 0.3)
#
# axes.set_aspect(1)
# axes.add_artist(draw_circle)
# plt.title('Circle')
# plt.show()
import numpy as np
import matplotlib.pyplot as plt

data = np.load('data_sub_diag_Bo7.npy')

x=data[0]
y=data[1]

fig=plt.figure()

#subplot1=fig.add_subplot(2,1,1)
plt.scatter(x, y)
plt.annotate('max value', xy=(100, 4), xytext=(2, 3),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )
# figure, axes = plt.subplots()
# draw_circle = plt.Circle((0.5, 0.5), 0.3,fill=False)
#
# plt.gcf().gca().add_artist(draw_circle)
# plt.title('Circle')
# axes.set_aspect(1)
plt.ylim(-1,0)
plt.show()
# subplot2=fig.add_subplot(2,1,2)
# subplot2.text(0.3, 0.5, '2nd Subplot')
#
# fig.suptitle("Add subplots to a figure")
