import numpy as np
import matplotlib.pyplot as plt
import math

def E(z):
    cg0 = -8 ** (0.2e1 / 0.3e1) * (math.pi * z) ** (0.2e1 / 0.3e1) * 4 ** (0.1e1 / 0.3e1) * (0.1e0 * z ** 3 + 5000 * z ** 2 + 0.125e2) ** (-0.2e1 / 0.3e1) / 8
    return cg0


def R(z,eta):
    # The exact fn R(z,eta) for example 1.
    # running the zero set through this function provides R(z,t=0)
    b = (1-np.cos(eta))
    a = -2*E(z)
    return (z/a)*b
    
zs = np.arange(0,5,0.1)





# data = np.genfromtxt("R_init.dat")

# This short block creates a file with the initial conditions for R. -> (z,R(z,t=0))
# outf = open("R_init.dat","w")
# for i in range(len(data[:,1])):
#     print(data[i,0],R(data[i,0],data[i,1]),file=outf)




# plt.plot(data[:,0],R(data[:,0],data[:,1]))
# plt.scatter(data[:,0],data[:,1])
# plt.show()

"""
plt.grid(True)
plt.xlim(0,5)
plt.ylim(0,6.28)
plt.plot(data[:,0],data[:,1])
plt.title("Initial conditions for R(z,t)")
plt.xlabel("z")
plt.ylabel("eta")
plt.show()

"""