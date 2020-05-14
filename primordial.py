import numpy as np
import math
"""
 Quasispherical Szekeres metric functions and arbitrary functions
 of integration defined for each subproblem

 Author:
 Nick Layden

"""
###############################################################################
# Quasispherical Szekeres functions
# Defined generally for the Szekeres metric

#################
## Dipole functions: S,P,Q
# including their z derivatives

def S(z):
    return z;

def Sp(z):
    return 1;

def P(z):
    return 0;

def Pp(z):
    return 0;

def Q(z): 
    return 0;

def Qp(z):
    return 0;

#################
## Dipole function: H
# including derivative H'

def H(z,x,y): 
	# Defined generally in terms of the dipole terms:
	# S, P, Q
    s = S(z);
    a = (x-P(z))/s;
    b = (y-Q(z))/s;
    h = s/2;
    return h*(1. + math.pow(a,2) + math.pow(b,2));


def Hprime(z,x,y):
	# H' defined generally in terms of the dipole functions:
	# S, P, Q, S', P', Q'

    a = x-P(z)
    b = y-Q(z)
    c = S(z)**2
    d = S(z)**3
    sp = Sp(z)
    pp = Pp(z)
    qp = Qp(z)

    e = 0.5*sp*(1 + a*a/c + b*b/c) + \
    	0.5*S(z)*(-2*a*pp/c - 2*a*a*sp/d - 2*b*qp/c - 2*b*b*sp/d )
    return e



###############################################################################
# Effective mass function
# And M_tilde which is scaled by H
# including derivative M'
# These terms differ by problem definitions
# Primordial black hole formation
def M(z):
	return (z**3)/2

def Mprime(z):
	return 3*(z**2)/2

def Mbar(z,x,y):
	# M_tilde function scaled by H
	return M(z)/H(z,x,y)

def Mbarprime(z,x,y):
	# Mb' = M' - 3MH'/H
    hp = Hprime(z,x,y)
    h = H(z,x,y)
    m = M(z)
    mp = Mprime(z)
    return mp - 3*m*hp/h

###############################################################################
# Energy function E(z) 
# Model D from Harada and Jhingan 2015
# Parameters n1,n2 > 2, rw=1, rc = 0.8-10
#
#
def E(z):
    return 0;



###############################################################################
# Functions subject to the field equations: Y, Y', density
# defined in terms of the field equation solutions: i.e. R, R', Rdot
# These generalized will depend on the dipole functions defined above
# and the choice of mass function.
#
# The functions will be defined in terms of t=const slices.
# This means that Y  = Y(R,z,x,y) 
#				  Y' = Y'(R,R',z,x,y)
#				 density = density(R,R',z.x,y)

def Y(R,z,x,y):
	# Y = R/H
	return R/H(z,x,y)

def Yprime(R,Rp,z,x,y):
	# Y' = (R' - RH'/H)/H
	a = 1/H(z,x,y)
	b = a*Hprime(z,x,y)
	c = a*(Rp - R*b)
	return c

def Density(R,Rp,z,x,y):
	# subject to the field equations
	# rho = 2M'/(8piY^2 Y')
	# where M is the 'scaled' Mass (Mbarprime)
	a = 2*Mbarprime(z,x,y)
	b = 8*np.pi
	c = math.pow(Y(R,z,x,y),2)
	d = Yprime(R,Rp,z,x,y)
	return a/(b*c*d)






























