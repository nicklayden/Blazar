import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("alantheta.dat")
data2 = np.genfromtxt("alan_f.dat")
# plt.plot(data[1,:],data[0,:],label="Theta func")
plt.plot(data[:,1],data[:,0], label="A(r)")
# plt.plot(data[:,3],data[:,2],c='r')
# plt.plot(data[:,5],data[:,4],c='r')
plt.grid(True)
plt.ylabel("r")
plt.xlabel("k")
plt.title("Zero set of $\sqrt{\Theta} = f(r)_{,r} + G(r)_{,r}	$")
plt.legend()
plt.show()