from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sb
import pandas as pd
import math as m
from primordial import *

fig = plt.figure()
ax = fig.add_subplot(111)
# ax = fig.add_subplot(111,projection='3d')

# data = np.genfromtxt("test.dat")
# R_init = np.genfromtxt("R_init.dat")
# t_vals = np.arange(0,18,0.01)
# single_curve = []
# all_curves = []


# def M(z):
#     return (z**3)/2.

# def Mprime(z):
#     return 1.5*z**2



def zero_2d(dat,x,y):
    for i in range(len(x)):
        for j in range(len(y)-1):
            aa = dat[i][j]
            bb = dat[i][j+1]
            if aa*bb <= 0:
                print("root!",i,j)
                ax.scatter(x_grid[i],y_grid[j],0)

def Rdot_e1(r,z,l):
    # Rdot for example 1. Check mass function.
    # let z = M
    return 

def R1_debnath(z,lam):

    a = 2./np.sqrt(lam);
    c = m.acos(-(3./2.)*np.sqrt(lam)*(2*M(z)));
    b = np.cos(c/3.);
    return a*b;

def R2_debnath(z,lam):
    a = 1./np.sqrt(lam);
    theta = m.acos(-(3./2.)*np.sqrt(lam)*(2*M(z)));
    b = -np.cos(theta/3.) + np.sqrt(3.) * np.sin(theta/3.);
    return a*b;


def horizon_cutoff(lam):
    # Determines z value of where horizons merge into 1 in this problem (example 2).
    return np.cbrt((2./3.)/np.sqrt(lam))

def horizon_cutoff_lam(z):
    return ((2./3.)/z**3)**2

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

def mu(rdot,r,z):
    e = E(z)
    a = np.sqrt(1 + 2*e);
    b = rdot + a;
    c = b/r;
    return -c/np.sqrt(2.0);

def rho(rdot,r,z):
    e = E(z)
    a = np.sqrt(1 + 2*e);
    b = rdot - a;
    c = b/r;
    return c/np.sqrt(2.0);

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

# time = np.arange(0,200,0.001)
time = np.genfromtxt("tsol.dat")
init = np.genfromtxt("z_out.dat")
mu_roots = np.genfromtxt("mu_zeros.dat")
app_roots = np.genfromtxt("apparent_zeros.dat")
Rprime_roots = np.genfromtxt("Rprime_zeros.dat")
R_init2d = np.genfromtxt("Example1/R_init.dat")
# init = np.linspace(1e-4, 1,50 )
for i in range(484):
    print(R_init2d[i][0], R_init2d[i][1])

surf = open("test2.dat",'r')
rprime = open("Rprime.dat",'r')
yprime = open("Yprime.dat")
rdot = open("rdot.dat",'r')
surf_d = []
rdot_d = []
yprime_d = []
rprime_d = []
for line in surf:
    # line.split()
    surf_d.append([line])
for line in rdot:
    rdot_d.append([line])
for line in rprime:
    rprime_d.append([line])
for line in yprime:
    yprime_d.append([line])


horz_r = []
horz_t = []
# fig = plt.subplots(111)

inv_r = []
inv_r2 = []
inv_t = []

# lam = m.pow(2./(3.*m.pow(1.0,3)),2);
lam = 1.0
# print("Maximum horizon formation point ",horizon_cutoff(lam))p_mu[:,dim2],c='r',lw=2.0,label="$\mu = 0 $")
# plt.plot(app_a[:,dim1],app)

split = 3

# identify the different horizons as a split between them
cosmo_mu = mu_roots[mu_roots[:,0] > split]
app_mu = mu_roots[mu_roots[:,0] < split]

# print(cosmo_mu[-1,1])

cosmo_a = app_roots[app_roots[:,0] > split]
app_a = app_roots[app_roots[:,0] < split]

# R = 0, z = 1, t = 2, t[i+1]=3, t[N]=4, R[0]=5,  t_c = 6
dim1 = 2
dim2 = 0

zmax = horizon_cutoff(lam)
print(zmax)
print(horizon_cutoff_lam(1.0))
# zs= np.linspace(0.7,0.8,100)
# R2h = [R2_debnath(i,lam) for i in zs]
# R1h = [R1_debnath(i,lam) for i in zs]
# plt.plot(zs,R2h,lw=3.0,c='b',label="$R_{CH}$")
# plt.plot(zs,R1h,lw=3.0,c='g',label="$R_{AH}$")


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.set_ylim3d(0,1)
R_initial = []
t_collapse = []
mu_surf = []
xg = 1
yg = 1

# for i in range(len(R_init2d)):
#     Ez = E(R_init2d[i][0])
#     raw = np.array(surf_d[i][0].split(),dtype=float)
#     rawp = np.array(surf_d[i][0].split(),dtype=float)
#     # mu_surf = []
#     # rho_surf = []
#     # for j in range(len(raw)):
#     #     mu_surf.append(mu(raw[j],rawp[j],Ez))
#     #     rho_surf.append(rho(raw[j],rawp[j],Ez))
#     # print(len(mu_surf))
#     # for j in range(len(raw)):
#     #     yprime_z.append(rawp[j] - raw[j]*Hprime(init[i],xg,yg)/H(init[i],xg,yg))
#     # scaled = raw
#     slice = np.array(raw,dtype=float)
#     # slice2= np.array(rho_surf,dtype=float)
#     # slice = np.array(surf_d[i][0].split(),dtype=float)

#     # t_collapse.append(time[len(slice)-2])
#     # R_initial.append(slice[0])
#     # ax.plot(time[0:len(slice)],slice,init[i],lw=2)
#     try:
#         plt.plot(time[0:len(slice)],slice,lw=0.5,c='b')#c=cm.gnuplot(init[i]))
#     except:
#         continue
#     # plt.plot(time[0:len(slice2)],slice2,lw=0.5,c='r')#c=cm.gnuplot(init[i])))
#     # try:
#     #     plt.plot(time[0:len(slice)],slice,lw=0.5,c='b')#c=cm.gnuplot(init[i]))
#     # #     # ax.plot(time[0:len(slice)],slice,init[i],lw=2)
#     # except:
#     #     continue



# Pick the solution along the curve z = 0.625... this is i = 50 in surf
# qq = 50
# t_index = -100

 
# shell_R  = np.array(surf_d[qq][0].split(),dtype=float)
# shell_Rp = np.array(rprime_d[qq][0].split(),dtype=float)
# shell_t  = time[0:len(shell_R)]
# print("z, t(coll - 10 steps), len(R_shell)    ")
# print(init[qq], shell_t[t_index], len(shell_R))

# zp = init[qq]
# Rr = shell_R[qq]
# Rpr = shell_Rp[qq]


# sc_radius = np.sqrt(m.pow(zp,2)*(Rr - zp*Rpr)/(zp*Rpr + Rr))
# print(sc_radius)

# l = 5*sc_radius
# ng = 100
# z_pos = init[qq]
# x_grid = np.linspace(-l,l,ng)
# y_grid = np.linspace(-l,l,ng)

# X,Y = np.meshgrid(x_grid,y_grid)

# kappa = 8*np.pi
# height = np.zeros((len(x_grid),len(y_grid)))
# ysol = np.zeros((len(x_grid),len(y_grid)))
# density = np.zeros((len(x_grid),len(y_grid)))
# for i in range(len(x_grid)):
#     for j in range(len(y_grid)):
#         xp = x_grid[i]
#         yp = y_grid[j]
#         zp = z_pos
#         R = shell_R[t_index]
#         Rp = shell_Rp[t_index]
#         height[i][j] = np.abs(Yprime(R,Rp,zp,xp,yp))
#         density[i][j] = Density(R,Rp,zp,xp,yp)

        # y_value = shell_R[t_index]/H(z_pos,x_grid[i],y_grid[j])

        # height[i][j] = (shell_Rp[t_index] - shell_R[t_index] * Hprime(z_pos,x_grid[i],y_grid[j])/H(z_pos,x_grid[i],y_grid[j]))/H(z_pos,x_grid[i],y_grid[j])
        # density[i][j]= m.pow(2*Mbarprime(z_pos,x_grid[i],y_grid[j])/(kappa*y_value**2 * height[i][j]),-1)
# rho = 2M'/(8pi Y^2 Y')




# df = pd.DataFrame(height,x_grid,y_grid)

# # fig = plt.figure()
# # ax = fig.add_subplot(111,projection='3d')
# sb.heatmap(df,vmin=-0.1,vmax=0.1,cmap='coolwarm',fmt='.1f',xticklabels=df.columns.values.round(2),
#                  yticklabels=df.index.values.round(2))
# ax.plot([0,1],[0,1])
# ax.set_xticks([])
# ax.set_yticks([])
# sb.set_style("ticks",{"xtick.major.size":8,"ytick.major.size":8})
# ax.plot_wireframe(X,Y,density)
# plt.colorbar()
# plt.plot(shell_t,shell_R)

# zero_2d(height,x_grid,y_grid)


# plt.plot(app_mu[:,1],app_mu[:,6])
# plt.plot(init,t_collapse,c='g',lw=2,label="Time until collapse")
# plt.plot(cosmo_mu[:,dim1],cosmo_mu[:,dim2],c='r',lw=2, label="$R_{CH}$ Formation")
# plt.plot(app_mu[:,dim1],app_mu[:,dim2],c='b',lw=2, label="$R_{AH}$ Formation")

# Apparent Horizon Detectors
plt.plot(app_mu[:,dim1],app_mu[:,dim2],c='r',lw=2.0,label="$\mu = 0 $")
plt.plot(app_a[:,dim1],app_a[:,dim2], ls=":",c='k',lw=2.0,label="$\Lambda R^3 + 6M - 3R = 0$")

# # # Cosmological Horizon Detectors
plt.plot(cosmo_mu[:,dim1],cosmo_mu[:,dim2],c='r',lw=2.0)# ,  label="$\mu = 0 $")
plt.plot(cosmo_a[:,dim1],cosmo_a[:,dim2], ls=":",c='k',lw=2.0)#, label="$\Lambda R^3 + 6M - 3R = 0$")

# # Shell Crossing Detector
# plt.plot(Rprime_roots[:,dim1],Rprime_roots[:,dim2],lw=2,c='r',label="$C_0 = 0$")




#     if len(slice) > 1:
#         for j in range(len(slice) - 1):
#             lhs = slice[j] - 2*init[i,0] 
#             rhs = slice[j+1] - 2*init[i,0] 
#             if (lhs * rhs) < 0 :
#                 # horz.append(time[j],slice[j])
#                 plt.scatter(time[j],slice[j],c='k')
#     for j in range(len(slice)):
#         slice[j] -= 2.0*init(j)
#     print(slice_)
#     print(time[0:len(slice_)])

# for i in range(len(surf_d)):
#     slice = np.array(surf_d[i][0].split(),dtype=float)
#     rdotslice = np.array(rdot_d[i][0].split(),dtype=float)
#     if len(slice) > 1:
#         for j in range(len(slice) - 1):
#             lhs = lambda_*(slice[j]**3) + 6*M(init[i]) - 3*slice[j]
#             rhs = lambda_*(slice[j+1]**3) + 6*M(init[i]) - 3*slice[j+1]

#             ml = rdotslice[j] 
#             mr = rdotslice[j+1]
#             if (ml * mr) < 0.:
#                 inv_r.append(slice[j])
#                 # inv_r.append(init[i])
#                 inv_t.append(time[j])
#                 continue

            # if (lhs * rhs) < 0 :
            #     # if init[i] > 0.11:
            #     horz_r.append(slice[j])
            #     # horz_r.append(init[i])
            #     horz_t.append(time[j])
            #     # plt.scatter(time[j],slice[j],c='k',marker=',')

# plt.plot(init,2*M(init))

# def ineq(lam):
#     return (2./3.)/np.sqrt(lam)
# lam = [ineq(lambda_) for i in range(len(init))]
# plt.plot(init,lam)

# print(len(horz_r))
# plt.plot(horz_t,horz_r,lw=3.0,c='g',label="Apparent Horizon")
# plt.scatter(inv_t,inv_r,c='r', label="Geometric Horizon")
# plt.scatter(inv_t,inv_r2,c='k')

# plt.scatter(mu_roots[:,2],mu_roots[:,0],c='g',lw=2.0,marker='+',label="$\mu = 0 $")
# plt.scatter(app_roots[:,2],app_roots[:,0],c='r',lw=1.0,marker='x',label="$\Lambda R^3 + 6M - 3R = 0$")


# plt.grid(True)
# plt.colorbar()
# plt.xlim(0,12.5)
# plt.ylim(0,2.8)
# plt.title(" Cosmological and Apparent Horizon Formation with $\Lambda = {0:2.3f} $".format(lam))

# plt.title("Collapse Time for a shell of Initial Radius R=R(z,t=0), $\Lambda = {0:2.3f} $".format(lam))

# plt.title("Map of $|Y_{,z}|$ with S(z) = $\sqrt{2}z$, $P(z)=Q(z)=4z $" + "\n $\Lambda = {0}$".format(lam))
plt.title("$\Lambda = 0.1$")


plt.grid(True)
plt.ylabel("R(t,z)")
plt.xlabel("time")
# plt.yscale("log")
plt.ylim(0,10)
plt.xlim(0,10)
plt.legend(loc="best")
# plt.colorbar()
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
