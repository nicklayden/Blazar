/**
 * Systems of equations for the Szekeres models
 * 
 * 
 * 
*/
#include <vector>
#include <math.h>

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


