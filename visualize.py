import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

def load_data(file):
	d = open(file,"r")
	data = []
	for line in d:
		row = []
		for item in line.split():
			row.append(float(item))
		data.append(row)
	return data

t = np.linspace(0,20,20001)

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


##########
# Construct dataframes and panels of the data
#
# convert None type values or NaN types to R=0 for R solutions
df_r = pd.DataFrame(data_rt,data_z,columns=t)
sb.plot_wireframe(df_r,robust=True)
plt.show()
# print(df_r)
 