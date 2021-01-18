import numpy as np
import matplotlib.pyplot as plt
import primordial as qs

		

test = qs.primordial_D(log=True)


z = np.arange(0.01,0.99,0.01)
y = [test.E(i) for i in z]
yp = [test.Eprime(i) for i in z]
ic = [test.init_cond(i) for i in z]
icp = [test.init_cond_prime(i) for i in z]

sc = [test.init_cond(i)/test.init_cond_prime(i) for i in z]

# plt.plot(z,y, label="E(z)")
# plt.plot(z,yp, label="E'(z)")
# plt.plot(z,ic, label="R(0,z)")
# plt.plot(z,icp, label="R'(0,z)")
plt.plot(z,sc, label="sc check")
plt.plot(z,z, label="pls below this line")
# plt.ylim(0,100)
plt.legend()
plt.grid(True)
plt.show()