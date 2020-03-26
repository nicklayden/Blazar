import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
from primordial import *



def load_data(file):
    # Loading solution curves from Blazar,
    # Data doesn't have consistent column sizes
    # Convert to float so pandas can read after
    d = open(file,"r")
    data = []
    for line in d:
        row = []
        for item in line.split():
            row.append(float(item))
        data.append(row)
    return data


def surface(f,xg,yg,R,Rp,z):
    # Compute the surface on the x-y plane 
    # for the function f
    # f(R,Rp,z,x,y) 
    result = np.zeros((len(xg),len(yg)))
    for i in range(len(xg)):
        for j in range(len(yg)):
            result[i][j] = f(R,Rp,z,xg[i],yg[j])
    return result


##########
# Loading all data from simulations, R', R, Rdot, Y', z, t
# Time grid currently is equivalent for all solution curves.
#   -- might need to change in the future.
data_rz = load_data("Rprime.dat")
data_yz = load_data("Yprime.dat")
data_rt = load_data("rdot.dat")
data_r  = load_data("test2.dat")
time    = np.genfromtxt("tsol.dat")
data_z  = np.genfromtxt("z_out.dat")

t = np.linspace(0,20,20001)


##########
# Construct dataframes and panels of the data
#
# convert None type values or NaN types to R=0 for R solutions
df_rt = pd.DataFrame(data_rt,data_z,columns=t)
df_rz = pd.DataFrame(data_rz,data_z,columns=t)
df_r  = pd.DataFrame(data_r ,data_z,columns=t)
 

#########
# Animation of the Y(t,r=c,x,y) function as a heatmap
#
# Indices for the z=constant surface and a starting point in time
t_i = 10
z_i = 50
a = 1
N = 100

zp = data_z[z_i]
R  = df_r.iloc[z_i][t[t_i]]
Rp = df_rz.iloc[z_i][t[t_i]]

x_g = np.linspace(-a,a,N)
y_g = np.linspace(-a,a,N)

y_surf = surface(Yprime, x_g,y_g,R,Rp,zp)

df_surf = pd.DataFrame(y_surf,y_g,columns=x_g)

sb.heatmap(df_surf)
plt.show()

# .iloc[i] grabs the ith row in a DataFrame
# for i in range(80):
#     plt.plot(t,df_r.iloc[i])
# plt.ylim(0,10)
# plt.show()


#########
# Plots of the data as heatmaps
#

# sb.heatmap(df_r.iloc[0:20],robust=True)
# plt.show()

# sb.heatmap(df_rt,robust=True)
# plt.show()













