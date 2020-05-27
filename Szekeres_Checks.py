# -*- coding: utf-8 -*-
"""
Created on Mon May 11 22:34:01 2020

@author: nicho
"""


import numpy as np
import matplotlib.pyplot as plt



def sc_check1(z):
    cg = 100 * z / (-2 * z ** 10 + z ** 8 + 1) ** 4 / (100 / (-2 * z ** 10 + z ** 8 + 1) ** 4 - 400 * z / (-2 * z ** 10 + z ** 8 + 1) ** 5 * (-20 * z ** 9 + 8 * z ** 7))
    return cg

def sc_check2(z):
    cg = 100 * z / (z ** 8 - 2 * z ** 6 + 1) ** 4 / (100 / (z ** 8 - 2 * z ** 6 + 1) ** 4 - 400 * z / (z ** 8 - 2 * z ** 6 + 1) ** 5 * (8 * z ** 7 - 12 * z ** 5))
    return cg

def Rprime(z):
    cg1 = -400 / (z ** 8 - 2 * z ** 6 + 1) ** 5 * (8 * z ** 7 - 12 * z ** 5)
    return cg1

z = np.arange(0,1,0.001)



plt.grid(True)

plt.xlim(0,1)
plt.ylim(0,1)

plt.xlabel("z")

plt.title("Initial Conditions with \n $r_c=10, r_w=1, n_1=8, n_2=6$")

plt.plot(z,z,label=r"$\frac{S}{S\prime}$")
plt.plot(z,sc_check2(z),label=r"$\frac{R}{R\prime}$")
plt.legend(loc='upper left')
plt.show()


plt.ylim(0,1)

plt.plot(z,Rprime(z))

plt.show()