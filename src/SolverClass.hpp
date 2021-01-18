/**
 *  ODE Solver class for use with Blazar program.
 *      
 *  Inputs: - System of ODEs. (As a class object) -- Using the boost::odeint framework
 *          - Initial Conditions. (Vector type)  
 *          - Time domain. (Vector type)  
 * 
 *  by: Nick Layden, 7 Dec 2020.
 * */
#include <vector>


class ODEBase {

    public: 
        ODEBase() {};
        // Pure virtual operator function. Must be defined in any inheriting classes.
        virtual void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) = 0;

};



class ODESolver: public ODEBase {

    public:
        ODESolver() {};



        // Input Variables
        std::vector<double> time_domain, initial_conditions;





};

















