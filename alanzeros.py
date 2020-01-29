import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("alantheta.dat")
data2 = np.genfromtxt("alan_f.dat")
# plt.plot(data[1,:],data[0,:],label="Theta func")
plt.plot(data2[1,:],data2[0,:], label="A(r)")

plt.grid(True)
plt.ylabel("r")
plt.xlabel("k")
plt.title("Zero set of $\sqrt{\Theta} = f(r)_{,r} + G(r)_{,r}	$")
plt.legend()
plt.show()