import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("alantheta.dat")
data2 = np.genfromtxt("alan_f.dat")
data3 = np.genfromtxt("para.dat")
plt.scatter(data3[:,1],data3[:,0],marker='.')
# plt.scatter(data3[:,3],data3[:,2],marker='.',c='r')
plt.scatter(data3[:,5],data3[:,4],marker='.',c='m')
for i in range(len(data[:,1])-1):
	print(data3[i,1]**2 + data3[i,0]**2)
# plt.plot(data[1,:],data[0,:],label="Theta func")
# plt.plot(data[:,1],data[:,0], label="A(r)")
# plt.plot(data[:,3],data[:,2],c='r')
# plt.plot(data[:,5],data[:,4],c='r')
plt.grid(True)
plt.ylabel("r")
plt.xlabel("k")
plt.title("Zero set of $\sqrt{\Theta} = f(r)_{,r} + G(r)_{,r}	$")
plt.legend()
plt.show()