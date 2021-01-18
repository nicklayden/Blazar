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

def E_model_D(z):
    # Parameters chosen to work for a model with Lambda=0,
    # Shell crossings appear INSIDE horizon.
    rc = 10.0
    rw = 1.0
    n1 = 8
    n2 = 9

    a = -0.5*math.pow(z/rc,2)
    b = math.pow(z/rw,n1)
    c = math.pow(z/rw,n2)
    d = math.pow(1 + b - 2*c,4)

    return a*d

# def Eprime(z):
#     cg = -z * math.pow(-0.2e1 * math.pow(z, 0.10e2) + math.pow(z, 0.8e1) + 0.1e1, 0.4e1) / 0.100e3 - z * z * math.pow(-0.2e1 * math.pow(z, 0.10e2) + math.pow(z, 0.8e1) + 0.1e1, 0.3e1) * (-0.20e2 * math.pow(z, 0.9e1) + 0.8e1 * math.pow(z, 0.7e1)) / 0.50e2;
#     return cg

def Eprime(z):
    # expanded version of Eprime above.
    rc = 10.0
    rw = 1.0
    n1 = 8
    n2 = 9

    a = -z/(rc**2)
    b = math.pow(1 + (z/rw)**n1 - 2*(z/rw)**n2,4)
    c = -2*(z/rc)**2
    d = math.pow(1 + (z/rw)**n1 - 2*(z/rw)**n2,3)
    e = n1*z**(n1-1)/(rw**n1) - 2*n2*z**(n2-1)/(rw**n2)

    return a*b + c*d*e

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


###############################################################################
# Cartan Invariants used as horizon detectors
#
# rho and mu
#
#

def mu(rd,r,E):
    # Cartan invariant mu. Collapsing phase horizon detector.
    # rd = R,t
    # r  = R
    # E  = E(z)
    a = rd - np.sqrt(1-2*E)
    return -np.sqrt(2)*a/(2*r)

################################################################################
#
#   Initial condition functions 
#   
#   For the functions R and R'
#

def e2_R0_max(z):
    # Based on the lambda=0 problem. maximum expansion 
    return -M(z)/E_model_D(z)

def example2_Rprime_init(z):
    Mp = 3*z*z/2.
    a = -Mp/E_model_D(z)
    b = (M(z)/math.pow(E_model_D(z),2))*Eprime(z)
    return a + b



class primordial_D():
    def __init__(self, log=False,n1=8,n2=9):
        # Initialization, need to set parameters for Model D for the E function
        self.n1 = n1
        self.n2 = n2
        self.rc = 10.0
        self.rw = 1.0

        if log:
            self.output()

    def E(self,z):
        # 'Energy' function E(z). z in (0,1). Defind by Model D in Harada/Jhingan
        a = -0.5*math.pow(z/self.rc,2)
        b = math.pow(z/self.rw,self.n1)
        c = math.pow(z/self.rw,self.n2)
        d = math.pow(1 + b - 2*c,4)

        return a*d

    def Eprime(self, z):
        # Derivative of E(z) function.

        n1 = self.n1
        n2 = self.n2
        rc = self.rc
        rw = self.rw

        a = -z/(rc**2)
        b = math.pow(1 + (z/rw)**n1 - 2*(z/rw)**n2,4)
        c = -2*(z/rc)**2
        d = math.pow(1 + (z/rw)**n1 - 2*(z/rw)**n2,3)
        e = n1*z**(n1-1)/(rw**n1) - 2*n2*z**(n2-1)/(rw**n2)

        return a*b + c*d*e


    def output(self):
        print("Model D parameters:")
        print("rc = ", self.rc)
        print("rw = ", self.rw)
        print("n1 = ", self.n1)
        print("n2 = ", self.n2)
        return

    def mass(self, z):
        # mass function M(z)
        return (z**3)/2.


    def mass_prime(self, z):
        # derivative a off M(z), M'(z)
        return 3*(z**2)/2.

    def init_cond(self, z):
        # Initial conditions for the problem, chosen to be 
        # R(0,z) = -M(z)/E(z)
        return -self.mass(z)/self.E(z)

    def init_cond_prime(self, z):
        # initial conditions for R'(0,z) (derivative of R(0,z)).
        m = self.mass(z)
        mp = self.mass_prime(z)
        e = self.E(z)
        ep = self.Eprime(z)

        return m*ep/(e**2) - mp/e


