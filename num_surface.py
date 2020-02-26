import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
import math as m
# data = np.genfromtxt("test.dat")
# R_init = np.genfromtxt("R_init.dat")
# t_vals = np.arange(0,18,0.01)
# single_curve = []
# all_curves = []
def M(z):
    return (z**3)/2.

def E(z):
    rc = 10.0;
    rw = 1.0;
    n1 = 8;
    n2 = 10;


    a = -0.5*m.pow(z/rc,2);
    b = m.pow(z/rw,n1);
    c = m.pow(z/rw,n2);
    d = m.pow(1 + b - 2*c,4);

    return a*d;

def Rdot(r,z,l):
    a = 2*E(z) + 2*M(z)/r + l*r*r/3.
    return -np.sqrt(a)

def mu(rdot,r,z,l):
    e = E(z)
    # rdot = Rdot(r,z,l)
    a = sqrt(1 + 2*e);
    b = rdot + a;
    c = b/r;
    return -c/np.sqrt(2.0);

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

time = np.arange(0,5000,0.001)
init = np.genfromtxt("z_out.dat")
# init = np.linspace(1e-4, 1,50 )

surf = open("test2.dat",'r')
# rdot = open("rdot.dat",'r')
surf_d = []
rdot_d = []
for line in surf:
    # line.split()
    surf_d.append([line])
# for line in rdot:
#     rdot_d.append([line])

horz_r = []
horz_t = []
# fig = plt.subplots(111)



lambda_ = 0.0
for i in range(len(surf_d)):
    slice = np.array(surf_d[i][0].split(),dtype=float)
    # rdotslice = np.array(rdot_d[i][0].split(),dtype=float)
    # print(len(slice))
    # print(len(time[0:len(slice)]))
    try:
        plt.plot(time[0:len(slice)],slice,lw=0.5,c='b')
        # plt.plot(time[0:len(slice)],mu(slice,init[i],lambda_))
        # plt.plot(time[0:len(rdotslice)],rdotslice,lw=0.5,c='r')
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

# for i in range(len(surf_d)):
#     slice = np.array(surf_d[i][0].split(),dtype=float)
#     # rdotslice = np.array(rdot_d[i][0].split(),dtype=float)
#     if len(slice) > 1:
#         for j in range(len(slice) - 1):
#             lhs = lambda_*(slice[j]**3) - 6*M(init[i]) - 3*slice[j]
#             rhs = lambda_*(slice[j+1]**3) - 6*M(init[i]) - 3*slice[j+1]

            # ml = rdotslice[j] 
            # mr = rdotslice[j+1]
            # if (ml * mr) <= 0.:
            #     horz_r.append(slice[j])
            #     horz_t.append(time[j])
            #     continue

            # if (lhs * rhs) < 0 :
            #     if init[i] > 0.11:
            #         horz_r.append(slice[j])
            #         # horz_r.append(init[i])
            #         horz_t.append(time[j])
            #         # plt.scatter(time[j],slice[j],c='k',marker=',')


# plt.scatter(horz_t,horz_r,lw=3.0,c='g',label="Geometric Horizon")


plt.grid(True)
# plt.colorbar()
# plt.xlim(0,12.5)
# plt.ylim(0,2.8)
# plt.title("Collapsing Phase with $\Lambda = 1.0$")
plt.ylabel("z (M)")
plt.xlabel("time")
# plt.legend()
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
