"""
    Tool for making plots of the zero sets of functions in python.

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")

# z_grid = np.genfromtxt("../z_grid.dat")
# eta_grid = np.genfromtxt("../eta_grid.dat")
# brackets = np.genfromtxt("../brackets.dat")
app_zeros = np.genfromtxt("../e1_R.dat")
rho_contracting = np.genfromtxt("../rho_ex1.dat")
# mu_contracting = np.genfromtxt("../mu_ex1.dat")
# surface = np.genfromtxt("../surface.dat")
# psim_zeros = np.genfromtxt("../psim.dat")
# Sort the columns by z value for plotting lines.
app_zeros[:,1] = app_zeros[:,1][np.argsort(app_zeros[:,0])]
app_zeros[:,0] = np.sort(app_zeros[:,0])

rho_contracting[:,1] = rho_contracting[:,1][np.argsort(rho_contracting[:,0])]
rho_contracting[:,0] = np.sort(rho_contracting[:,0])

# mu_contracting[:,1] = mu_contracting[:,1][np.argsort(mu_contracting[:,0])]
# mu_contracting[:,0] = np.sort(mu_contracting[:,0])

# print(np.shape(surface))
# surface_masked = np.ma.masked_where( surface > -1 , surface )
# print(np.log10(surface))
# surface[np.isnan(surface)] = 0
# print(np.shape(surface))
# X,Y = np.meshgrid(z_grid,eta_grid)



def tC(x):
    return 0.1*x**3 + 12.5
def tB(x):
    return -5000.*x**2

xvals = np.arange(0,10,0.10)


plt.title(r"Zero set of $\mu $ on top of AH at R=2M")
plt.ylabel("t  ")
plt.xlabel("M ")
plt.plot(app_zeros[:,0],app_zeros[:,1],label="R=2M")
# plt.plot(mu_contracting[:,0],mu_contracting[:,1],label="$\mu$")
plt.plot(rho_contracting[:,0],rho_contracting[:,1],label=r" $ \mu  (contracting)$ ",color="r",ls="--")
# plt.plot(xvals,tC(xvals), label="tC(M)")
# plt.plot(psim_zeros[:,0],psim_zeros[:,1],label="$\psi_-$")
plt.ylim(10,25)
plt.xlim(0,7)
plt.grid(True)
plt.legend(loc="upper right")
plt.show()







# ax.plot_surface(X,Y,surface[80:100,0:99],alpha=0.5)
# ax.plot(zeros[:,0],zeros[:,1],0,color="r")
# ax.plot(brackets[:,0],brackets[:,1],0,color="g")
# ax.plot(brackets[:,2],brackets[:,3],0,color="g")
# ax.set_zlim(-1,1)
# plt.show()

# plt.imshow(surface)

# plt.show()






