import numpy as np
from scipy.integrate import quad
import math as m
# double R1_debnath(double z, double lambda) {
#     double a,b,c;
#     a = 2./sqrt(lambda);

#     c = acos(-(3./2.)*sqrt(lambda*2*mass_e2(z)));
#     b = cos(c/3.);
#     return a*b;
# }

# double R2_debnath(double z, double lambda) {
#     double a,b,c,theta;
#     a = 1./sqrt(lambda);
#     theta = acos(-(3./2.)*sqrt(lambda*2*mass_e2(z)));
#     b = -cos(theta/3.) + sqrt(3.) * sin(theta/3.);
#     return a*b;
# }

def E(z):
    rc = 10.0
    rw = 1.0
    n1 = 8
    n2 = 10


    a = -0.5*m.pow(z/rc,2)
    b = m.pow(z/rw,n1)
    c = m.pow(z/rw,n2)
    d = m.pow(1 + b - 2*c,4)

    return a*d;

def f(R,z,lam):
	g = m.pow(z,3)
	a = 2*E(z)
	b = g/R
	c = lam*m.pow(R,2)/3.

	return 1./np.sqrt(a + b + c)

def Rh1(z,lam):
	M = m.pow(z,3)/2.
	a = 2./np.sqrt(lam)
	c = m.acos(-(3./2.)*np.sqrt(lam*M))
	b = np.cos(c/3.)
	return a*b;


print(Rh1(0.62482899999999997,1))