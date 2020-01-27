import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("alanroots.dat")
plt.plot(data[1,:],data[0,:])

plt.grid(True)
plt.ylabel("r")
plt.xlabel("k")
plt.title("Zero set of $\sqrt{\Theta} = f(r)_{,r} + G(r)_{,r}	$")

plt.show()