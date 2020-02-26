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
        double m_z,m_f;
        ode_e1(double z, bool signflag): m_z(z),m_f(signflag) {};
        void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) {
            double R,Lambda;
            R = xin[0];
            Lambda=1.0;
            // z values are independent of each other. so the solution will be evolved for each z.
            // The value being integrated here is R.
            
            // Positive root
            if(m_f) { 
            dxdt[0] = -sqrt(2*E(m_z) + 2*m_z/R);
            } else { 
            // Negative root
            dxdt[0] = -sqrt(2*E(m_z) + 2*m_z/R + Lambda*R*R/3.);         
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
        ode_e2(double z, double lambda,  bool signflag): m_z(z),m_l(lambda),m_f(signflag) {};
        double m_z,m_l;
        bool m_f;
        std::vector<double> full_sol,mu_sol;
        // std::ofstream stats_out.open("stats.dat");
        void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) {
            double R,mu,rho;
            R = xin[0];
            // z values are independent of each other. so the solution will be evolved for each z.
            // The value being integrated here is R.
            
            // Positive root
            if(m_f) { 
                dxdt[0] = sqrt(2*E(m_z) + 2*M(m_z)/R  + m_l*R*R/3.);
            } else { 
                // Negative root
                dxdt[0] = -sqrt(2*E(m_z) + 2*M(m_z)/R + m_l*R*R/3.);
                // xin[1] = dxdt[0];
                // mu = mu_example2(dxdt[0],R,E(m_z));
                // rho = rho_example2(dxdt[0],R,E(m_z));
                // xin[2] = mu;
                // xin[3] = rho;
                // std::cout << mu << std::endl;
                // std::cout << dxdt[0] << " " << R << " " << mu << " " << rho << std::endl;
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
        }


};

