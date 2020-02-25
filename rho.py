import numpy as np
import matplotlib.pyplot as plt
import math as m

def energy(z):
    rc = 10.0
    rw = 1.0
    n1 = 8
    n2 = 10


    a = -0.5*m.pow(z/rc,2)
    b = m.pow(z/rw,n1)
    c = m.pow(z/rw,n2)
    d = m.pow(1 + b - 2*c,4)

    return a*d

def mass(z):
	return m.pow(z,3)/2.

def Rmax(z):
	return -mass(z)/energy(z)


z = 0.8001


data = np.genfromtxt("mu.dat")
# mud = np.genfromtxt("mu.dat")
rho = data[1::6,3]
mu = data[1::6,2]
rdot = data[1::6,0]
r = data[1::6,1]
# time2 = np.linspace(0,1000,len(mud))
time = np.arange(0,1000,0.01)

lam = 0.
geom = []
for i in range(len(r)):
	lhs = lam * (r[i]**3) - 6*mass(z) - 3*r[i]
	geom.append(lhs)
	# rhs = lam * (r[i+1]**3) - 6*mass(z) - 3*r[i+1]

	# if lhs*rhs <= 0.:
	# 	print(lhs, time[i])

	# lhsm = mu[i]
	# rhsm = mu[i+1]

	# if lhsm*rhsm <= 0.:
	# 	print(mu[i], time[i])



plt.plot(time,geom,label="geom")
print(len(time), len(rho))
plt.plot(time, rho, label="rho")
# plt.plot(time2, mud, label="other mu")
plt.plot(time,mu, label="mu")
plt.plot(time,r, label="r")
plt.plot(time,rdot,label="rdot")

plt.legend()
plt.grid(True)
plt.show()
