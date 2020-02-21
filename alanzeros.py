import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("alan_f.dat")
# data = np.genfromtxt("alan_f.dat")
data3 = np.genfromtxt("para.dat")
# data4 = np.genfromtxt("eta_bang_f.dat")
theta = np.arange(0,2*np.pi,0.01)
# plt.scatter(np.cos(theta),np.sin(theta))
# plt.scatter(data3[:,1],data3[:,0],marker='.')
# # plt.scatter(data3[:,3],data3[:,2],marker='.',c='r')
# plt.scatter(data3[:,5],data3[:,4],marker='.',c='m')
# for i in range(len(data[:,1])-1):
# 	print(data3[i,1]**2 + data3[i,0]**2)
data = data[data[:,0].argsort()]

plt.plot(data[:,1],data[:,0],c='r')
plt.scatter(data[:,3],data[:,2],marker='.')
plt.scatter(data[:,5],data[:,4],marker='.')
# plt.scatter(data4[:,1],data4[:,0], label="A(r)",marker='.')
# plt.scatter(data3[:,3],data3[:,2],c='r')
# plt.scatter(data3[:,5],data3[:,4],c='r')
plt.grid(True)
plt.ylabel("r")
plt.xlabel("k")

# plt.title("Zero set of $\sqrt{\Theta} = f(r)_{,r} + G(r)_{,r}	$")
# plt.legend()
plt.show()