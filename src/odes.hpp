/**
 * Systems of equations for the Szekeres models
 * 
 * 
 * 
*/
#include <vector>
#include <math.h>
#include <fstream>


struct push_back_state_and_time 
{
    std::vector<std::vector<double> >& m_states;
    std::vector<double>& m_times;

    push_back_state_and_time(std::vector<std::vector<double> >& states, std::vector<double>& times)
    :m_states(states), m_times(times) { }

    void operator()(const std::vector<double>& x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }


};

/** 
 * EXAMPLE 1: Galactic Black Hole Formation
 * 
 * Coordinates are changed to z=M. So the z coordinate here is actually the mass
 * ranges on values are z = [0,inf), R = [0,inf) 
 * 
 * Initial conditions retrieved from numerically solving the exact solution at t=0
*/

class ode_e1 {
    public: 
        double m_z,m_f, m_l;
        ode_e1(double z, double lam, bool signflag): m_z(z),m_l(lam),m_f(signflag) {};
        void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) {
            double R;
            R = xin[0];
            // z values are independent of each other. so the solution will be evolved for each z.
            // The value being integrated here is R.
            
            // Positive root
            if(m_f) { 
            dxdt[0] = -sqrt(2*E(m_z) + 2*m_z/R + m_l*R*R/3.);
            } else { 
            // Negative root
            dxdt[0] = -sqrt(2*E(m_z) + 2*m_z/R + m_l*R*R/3.);         
            }

        }

        double E(double z) {
            return -pow(0.8e1, 0.2e1 / 0.3e1) * pow(0.3141592654e1 * z, 0.2e1 / 0.3e1) * pow(0.4e1, 0.1e1 / 0.3e1) * pow(0.1e0 * pow(z, 0.3e1) + 0.5000e4 * z * z + 0.125e2, -0.2e1 / 0.3e1) / 0.8e1;
        }

};


/** 
 * EXAMPLE 2: Primordial Black Hole Formation
 *
 * Initial conditions retrieved from numerically solving the exact solution at t=0
 * Big bang time tB(z) = 0 for this problem.
*/

class ode_e2 {
    public: 
        ode_e2(double z, double lambda,  bool signflag, double x=0, double y=0): m_z(z),m_l(lambda),m_f(signflag),x(x),y(y) {};
        double m_z,m_l;
        // Fixing to the x=y=0 plane for now.
        double x,y;
        bool m_f;
        std::vector<double> full_sol,mu_sol;
        // std::ofstream stats_out.open("stats.dat");
        void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) {
            double R,mu,rho, Rp;
            R = xin[0];
            Rp = xin[1];
            // z values are independent of each other. so the solution will be evolved for each z.
            // The value being integrated here is R.
            
            // Positive root
            if(m_f) { 
                dxdt[0] = sqrt(2*E(m_z) + 2*M(m_z)/R  + m_l*R*R/3.);
                dxdt[1] = udot(Rp,dxdt[0],R,m_z,m_l);
            } else { 
                // Negative root: Rdot equation
                dxdt[0] = -sqrt(2*E(m_z) + 2*M(m_z)/R + m_l*R*R/3.);
                dxdt[1] = udot(Rp,dxdt[0],R,m_z,m_l);
                

                // WRONG
                // Ydot equations -- full quasispherical problem here.
                // dxdt[0] = -sqrt(2*M(m_z)/(R*pow(H(m_z,x,y),3)) + 2*E(m_z)/(pow(H(m_z,x,y),2)) + m_l*R*R/3. );
                // dxdt[1] = udot_y(Rp,dxdt[0],R);

                // If solving for Y now, need to rescale to get R back!

            }

            // full_sol.push_back(dxdt[0]);
        }

        double E(double z) {
            double a,b,c,d;
            double rc,rw;
            int n1, n2;
            rc = 10.0;
            rw = 1.0;
            n1 = 8;
            n2 = 10;


            a = -0.5*pow(z/rc,2);
            b = pow(z/rw,n1);
            c = pow(z/rw,n2);
            d = pow(1 + b - 2*c,4);

            return a*d;
            // return 0.0;
        }

        double rho_example2(double rdot, double r, double energy) {
            // Extended Cartan invariant that detects the horizon
            double a,b,c,d;
            a = sqrt(1 + 2*energy);
            b = rdot - a;
            c = b/r;
            return c/sqrt(2.0);
        }

        double mu_example2(double rdot, double r, double energy) {
            // Extended Cartan invariant that detects the horizon
            double a,b,c,d;
            a = sqrt(1 + 2*energy);
            b = rdot + a;
            c = b/r;
            return -c/sqrt(2.0);
        }

        double M(double z) {
            // M = z^3/2   Jhinghan paper
            return z*z*z/2.;

            // de Sitter Universe M(r) = 0;
            // return 0;
            // double lambda = 0.1;
            // return (1./3.)/sqrt(lambda);
        }

        double udot(double u, double rd, double r,double z, double lambda) {
            // Rprime dot equation, trying to solve for Rprime for each time in the domain.
            // This helps calculate rho the energy density.
            // This is a second ODE to solve concurrently.
            double a,b,c,d,e;
            a = 1./rd;
            b = Mprime(z)/r;
            c = -M(z)*u/pow(r,2.);
            d = Eprime(z);
            e = lambda*r*u/3.;
            return a*(b+c+d+e);


        }

        double udot_y(double u, double rd, double r) {
            // set fixed values
            double m,mp,h,hp,e,ep;
            double a,b,c,d,f,g;
            m = M(m_z);
            mp = Mprime(m_z);
            h = H(m_z,x,y);
            hp = Hprime(m_z,x,y);
            e = E(m_z);
            ep = Eprime(m_z);

            // function terms
            a = mp/(pow(h,3)*r);
            b = -3*m*hp/(pow(h,4)*r);
            c = -m*u/(pow(h,3)*pow(r,2));
            d = ep/pow(h,2);
            f = -2*e*hp/pow(h,3);
            g = m_l*r*u/3;

            return (a+b+c+d+f+g)/rd;   
        }

        double H(double z,double x, double y) {
            double h,a,b,s;
            s = S(z);
            a = (x-P(z))/s;
            b = (y-Q(z))/s;
            h = s/2;
            return h*(1. + pow(a,2) + pow(b,2));
        }
        
        double Hprime(double z, double x, double y) {
            return (Sprime(z)*(pow(S(z),2) - x*x - y*y))/(2*pow(S(z),2));
        }

        double Sprime(double z) {
            return sqrt(2);
        }

        double S(double z) {
            return sqrt(2)*z;
        }

        double P(double z) {
            return 0;
        }

        double Q(double z) {
            return 0;
        }



        double Mprime(double z) {
            return 3.*z*z/2.;
        }
        double Eprime(double z) {
            double cg;
            cg = -z * pow(-0.2e1 * pow(z, 0.10e2) + pow(z, 0.8e1) + 0.1e1, 0.4e1) / 0.100e3 - z * z * pow(-0.2e1 * pow(z, 0.10e2) + pow(z, 0.8e1) + 0.1e1, 0.3e1) * (-0.20e2 * pow(z, 0.9e1) + 0.8e1 * pow(z, 0.7e1)) / 0.50e2;
            return cg;
        }

};

