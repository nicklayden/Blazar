import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

# data = np.genfromtxt("test.dat")
# R_init = np.genfromtxt("R_init.dat")
# t_vals = np.arange(0,18,0.01)
# single_curve = []
# all_curves = []
def M(z):
    return (z**3)/2.

# data2 = []
# infile = open("test2.dat","r")
# for row in infile:
#     data2.append(row)

# print(data[1])

# for i in range(len(data[:,0])):
#     for j in range(len(data[0,:])):
#         # print(data[i,j], np.isnan(data[i,j]) ,end=' ')
#         if(~np.isnan(data[i,j])):
#             single_curve.append(data[i,j])
#         # else:
#             # single_curve.append(99)
#     all_curves.append(single_curve)

# print("Removed NaN Values.")

# ax = plt.axes(projection='3d')
# ax.plot3D(all_curves[100],np.arange(0,len(all_curves[100]),1))
# plt.show()

time = np.arange(0,1000,0.01)
init = np.genfromtxt("z_out.dat")
# init = np.linspace(1e-4, 1,50 )
print(len(init))
surf = open("test2.dat",'r')
surf_d = []
for line in surf:
    # line.split()
    surf_d.append([line])

print(len(surf_d))
horz_r = []
horz_t = []
# fig = plt.subplots(111)

for i in range(len(surf_d)):
    slice = np.array(surf_d[i][0].split(),dtype=float)
    # print(len(slice))
    try:
        plt.plot(time[0:len(slice)],slice,lw=0.5,c=cm.jet(init[i]))
    except:
        continue
    # if len(slice) > 1:
    #     for j in range(len(slice) - 1):
    #         lhs = slice[j] - 2*init[i,0] 
    #         rhs = slice[j+1] - 2*init[i,0] 
    #         if (lhs * rhs) < 0 :
    #             # horz.append(time[j],slice[j])
    #             plt.scatter(time[j],slice[j],c='k')
    # for j in range(len(slice)):
    #     slice[j] -= 2.0*init(j)
    # print(slice_)
    # print(time[0:len(slice_)])
lambda_ = 0.001
for i in range(len(surf_d)):
    slice = np.array(surf_d[i][0].split(),dtype=float)
    if len(slice) > 1:
        for j in range(len(slice) - 1):
            lhs = lambda_*(slice[j]**3) - 6*M(init[i]) - 3*slice[j]
            rhs = lambda_*(slice[j+1]**3) - 6*M(init[i]) - 3*slice[j+1]
            if (lhs * rhs) < 0 :
                if init[i] > 0.11:
                    # horz_r.append(slice[j])
                    horz_r.append(init[i])
                    horz_t.append(time[j])
                    # plt.scatter(time[j],slice[j],c='k',marker=',')


plt.plot(horz_t,horz_r,lw=3.0,c='r',label="Geometric Horizon")


plt.grid(True)
# plt.colorbar()
# plt.xlim(0,12.5)
# plt.ylim(0,2.8)
plt.title("Collapsing Phase with $\Lambda = 1.0$")
plt.ylabel("z (M)")
plt.xlabel("time")
plt.legend()
plt.show()
# print(np.array(len(surf_d[20][0].split(" "))))
# test = np.array(len(surf_d[40][0].split(" ")))


# exact = np.genfromtxt("output/e1_R.dat")
# data = np.genfromtxt("eta_bang.dat")
# plt.grid(True)

# newgrid = np.arange(0,5,0.01)
# yinterp = np.interp(newgrid, data[:,0],data[:,1])




# exact[:,1] = exact[:,1][np.argsort(exact[:,0])]
# exact[:,0] = np.sort(exact[:,0])



# plt.plot(data[:,0],data[:,1],linewidth=3.0,label="Numerical solution")
# plt.scatter(newgrid,yinterp,c="r")
# plt.scatter(data[:,0],data[:,1],linewidth=3.0,label="Numerical solution")
# plt.scatter(data[208:,2],data[208:,3])
# plt.plot(exact[:,0],exact[:,1],ls="--",linewidth=3.0, label="Exact solution")
# plt.xlim(0,4.5)
# plt.ylim(10,20)
# plt.xlabel("z")
# plt.ylabel("t")
# plt.legend()
# plt.show()
