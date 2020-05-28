/**
 *  Program for solving for the zero set of a function of 2 variables.
 * 
 *  Made for analyzing the zero sets of Cartan Invariants
 * 
 *  Project reconfigured into "Blazar" as a git repo
 *  author: Nick Layden
 * 
*/
#include <GL/glut.h>

// Standard libraries
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#include <complex>
#include <chrono>
#include <limits>
// #include <gmpxx.h>

// Boost libraries
#include <boost/multiprecision/mpfr.hpp>
// Note: in mpfr, expression templates are disabled. This reduces compiles but increases runtime.
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint.hpp>

// Graphics libraries
// #include <SFML/Graphics.hpp>


// Custom headers
#include "methods.hpp"
#include "invariants.hpp"
#include "odes.hpp"
#include "testfunctions.hpp"

// Type definitions
// typedef boost::multiprecision::cpp_dec_float_50 mp_type;
namespace mp = boost::multiprecision;
typedef mp::mpfr_float_100 mp_type; // faster compile time


// Preprocessor defines
#define PI M_PIl

void compute_and_save_zero_set(mp_type (*f)(mp_type,mp_type), std::vector<std::vector<mp_type> > domain, std::string filename);
std::vector<std::vector<double> > read_file(std::ifstream& input);



inline double S(double z) {
    return sqrt(2)*z;
}

inline double Sprime(double z) {
    return sqrt(2);
}

inline double P(double z) {
    return 0;
}

inline double Q(double z) {
    return 0;
}

inline double H(double z,double x, double y) {
    double h,a,b,s;
    s = S(z);
    a = (x-P(z))/s;
    b = (y-Q(z))/s;
    h = s/2;
    return h*(1. + pow(a,2) + pow(b,2));
}

inline double Hprime(double z, double x, double y) {
    return (Sprime(z)*(pow(S(z),2) - x*x - y*y))/(2*pow(S(z),2));
}

inline double Eprime(double z) {
    double cg;
    cg = -z * pow(-0.2e1 * pow(z, 0.10e2) + pow(z, 0.8e1) + 0.1e1, 0.4e1) / 0.100e3 - z * z * pow(-0.2e1 * pow(z, 0.10e2) + pow(z, 0.8e1) + 0.1e1, 0.3e1) * (-0.20e2 * pow(z, 0.9e1) + 0.8e1 * pow(z, 0.7e1)) / 0.50e2;
    return cg;
}

inline double Rprime(double z) {
    double cg0;
    cg0 = -0.100e3 * (0.78e2 * pow(z, 0.10e2) - 0.31e2 * pow(z, 0.8e1) + 0.1e1) * pow(0.2e1 * pow(z, 0.10e2) - pow(z, 0.8e1) - 0.1e1, -0.5e1);
    return cg0;
}


inline mp_type rho_example2(mp_type rdot, mp_type r, mp_type energy) {
    // Extended Cartan invariant that detects the horizon
    mp_type a,b,c,d;
    a = sqrt(1 + 2*energy);
    b = rdot - a;
    c = b/r;
    return c/sqrt(2.0);
}

inline mp_type alan_f(mp_type r, mp_type k) {
    // Alans specific prob. Just wanted a graph showing the function
    // crossing f(r,k)=0 for r near pi/4 and k = (1,2)

    // NOte: Function is likely written incorrectly in paper.
    mp_type a,b,c,d,e,f,g,h;
    a = sinh(k*(M_PIl - r)) + sinh(k*r);
    b = cosh(k*(M_PIl - r)) - cosh(k*r);
    c = sinh(k*(M_PIl/2. + r)) + sinh(k*(M_PIl/2. - r));
    d = cosh(k*(M_PIl/2. + r)) - cosh(k*(M_PIl/2. - r));

    e = -(1./sin(r));
    f = -cos(r)/(sin(r)*sin(r));
    g = 1./cos(r);
    h = tan(r);

    return e*(  f*a + k*b) + g*( h*c + k*d )  ;
}


inline mp_type alan_f_A(mp_type r, mp_type k) {
    // Alans function
    // THIS IS A(r)  -- equation 86/87 in paper
    mp_type a,b,c,d,e,f;
    a = sinh(k*M_PI);
    b = cosh(k*r);
    c = cosh(k*M_PI);
    d = sinh(k*r);
    e = -cos(r)/(sin(r)*sin(r));
    f = k/sin(r);

    return e*(a*b - c*d + d) + f*(a*d - c*b + b);
}

inline mp_type alan_f_Gr(mp_type r, mp_type k) {
    // Alan's function G(r),r
    // Depends on alan_f_A for calculating A(r2)
    mp_type a,b,c,d,r2;

    a = sqrt(2.)*cos(r);
    b = sqrt(1. - 0.5*sin(r)*sin(r));
    c = a/b;
    r2 = acos(sin(r)*(-1/sqrt(2.)));

    d = alan_f_A(r2,k);
    return c*d;
}

inline mp_type alan_const_k(mp_type r, mp_type k) {
    return alan_f_A(r,1.0);
}


inline mp_type alan_theta(mp_type r, mp_type k) {
    // Alan's function sqrt(theta) = F,r + G,r
    // Depends on alan_f_Gr, alan_f
    return alan_f_A(r,k) + alan_f_Gr(r,k);
}

inline mp_type alan_ar_test(mp_type r, mp_type k) {
    return alan_f_A(r,k) - alan_f_A(M_PI/2. - r, k); 
}

inline mp_type alan_ar_test2(mp_type r, mp_type k) {
    return alan_f_A(M_PI/2. - r, k); 
}

inline mp_type test_paraboloid(mp_type x, mp_type y) {
    // roots are the circle x^2 + y^2 = 1
    return x*x + y*y - 1.0;
}

inline mp_type plane(mp_type x, mp_type y) {
    //eqn for a plane: ax + by + cz = d
    // as f(x,y) -> z = (d - ax - by)/c
    mp_type a,b,c,d;
    a = 1.0;
    b = 10.0;
    c = 2.0;
    d = 1.0;

    return (d - a*x - b*y)/c;
}

inline double energy_e2(double z) {
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

inline double energy_e11(double z) {
    return -pow(0.8e1, 0.2e1 / 0.3e1) * pow(0.3141592654e1 * z, 0.2e1 / 0.3e1) * pow(0.4e1, 0.1e1 / 0.3e1) * pow(0.1e0 * pow(z, 0.3e1) + 0.5000e4 * z * z + 0.125e2, -0.2e1 / 0.3e1) / 0.8e1;
}

inline double mass_e2(double z) {
    return pow(z,3)/2.;
    // double lambda = 0.1;
    // return (1./3.)/sqrt(lambda);
}

inline double example2_eta_set(double z, double eta) {
    return eta - sin(eta) + 10*pow(-2*energy_e2(z),3./2.)/mass_e2(z) ;
}

inline double example2_Rmax(double z) {
    // Start the areal radii at their max value so we can evolve the contracting phase (negative root ODE)
    return -mass_e2(z)/energy_e2(z);
}

inline double example2_Rprime_init(double z) {
    double a,b, Mp;
    Mp = 3*z*z/2.;
    a = -Mp/energy_e2(z);
    b = (mass_e2(z)/pow(energy_e2(z),2))*Eprime(z);
    return a + b;
}


inline double Rdot(double r, double z, double lambda) {
    // CGECK THIS
    double a;
    a = 2*energy_e2(z) + 2*mass_e2(z)/r + lambda*r*r/3.;
    return -sqrt(a);
}


inline double Rdot_e1(double r, double z, double lambda) {
    // CGECK THIS
    double a;
    a = 2.*energy_e11(z) + (2.*z)/r + lambda*r*r/3.;
    return -sqrt(a);
}


double mu_example2(double rdot, double r, double energy) {
    // Extended Cartan invariant that detects the horizon.
    // Generalized to fit examples 1 and 2.
    double a,b,c;
    a = sqrt(1 + 2*energy);
    b = rdot + a;
    c = b/r;
    return -c/sqrt(2.0);
}

double rho_example2(double rdot, double r, double energy) {
    // Extended Cartan invariant that detects the horizon
    double a,b,c;
    a = sqrt(1 + 2*energy);
    b = rdot - a;
    c = b/r;
    return c/sqrt(2.0);
}

double apparent_horizon(double R, double z, double lambda) {
    // Debnath Nath Chakraborty 2006 paper equaiton 25
    // Requires mass for EXAMPLE 2
    return lambda*pow(R,3) + 6*mass_e2(z) - 3*R;
}


double apparent_horizon_e1(double R, double z, double lambda) {
    // Debnath Nath Chakraborty 2006 paper equaiton 25
    // Requires mass for EXAMPLE 2
    return lambda*pow(R,3) + 6*z - 3*R;
}


double integrand(double R, double z, double lambda) {
	// Integral of evolution equation for t_collapse time

	double a,b,c,d;

	a = 2*mass_e2(z)/R;
	b = 2*energy_e2(z);
	c = lambda*pow(R,2)/3.;
	d = a+b+c;

	return 1./sqrt(d);


}

double zero_force_init(double z, double lambda) {
	


	return pow(mass_e2(z)/(3.*lambda)  ,1./3.);
}





double trapz(double Rh, double z, double lambda) {
	// create domain in R
	int N = 10000;
	double dr = Rh/N;
	double sum = 0.0;
	std::vector<double> Rdomain;
	std::vector<double> answer;
	for (int i = 0; i < N; ++i)
	{
		Rdomain.push_back(i*(Rh)/N);
	}


	std::vector<double> integral_parts;
	for (int i = 0; i < N; ++i)
	{
		integral_parts.push_back(integrand(Rdomain[i],z,lambda));
		// integral_parts.push_back(R_domain[i]*Rdomain[i]);
	}

	for (int i = 1; i < N-1; ++i)
	{
		sum += 2*integral_parts[i];
	}

	sum += integral_parts[0];
	sum += integral_parts[N-1];

	return dr*sum/2.;
}






void zeros_output(std::vector<std::vector<double> > rsol, std::vector<std::vector<double> > function,std::vector<double> z, std::vector<double> t, std::string file, double lambda) {
    // Compute simple bracketing set , return one half as a pseudo zero set
    std::vector<std::vector<double> > zeros;
    std::vector<double> slice;
    double lhs,rhs;
    int jmax;

    for (int i = 0; i < function.size(); ++i)
    {
        slice.clear();
        for (int j = 0; j < function[i].size()-1; ++j)
        {
            jmax = function[i].size()-1;
        	// Compute a bracketing set for the zeros of the function in question
            lhs = function[i][j];
            rhs = function[i][j+1];
            if (lhs * rhs <= 0.0) {
                slice.push_back(rsol[i][j]); // left bracket of the zero 
                slice.push_back(z[i]); // z position of the zero
                slice.push_back(t[j]); // time of the zero
                slice.push_back(t[j+1]); // time bracket of the zero
                slice.push_back(t[jmax]); // maximum time until solution breaks 
                slice.push_back(rsol[i][0]); // initial value of R
                slice.push_back(trapz(rsol[i][0],z[i],lambda)); // collapse time integral from R_0

                zeros.push_back(slice);
                slice.clear();
            }
            

        }
    }
    matrix_to_file3(zeros,file);
}


double R1_debnath(double z, double lambda) {
    double a,b,c;
    a = 2./sqrt(lambda);

    c = acos(-(3./2.)*sqrt(lambda)*(2*mass_e2(z)));
    b = cos(c/3.);
    return a*b;
}

double R2_debnath(double z, double lambda) {
    double a,b,c,theta;
    a = 1./sqrt(lambda);
    theta = acos(-(3./2.)*sqrt(lambda)*(2*mass_e2(z)));
    b = -cos(theta/3.) + sqrt(3.) * sin(theta/3.);
    return a*b;
}

inline double E2_init_focus_r(double z, double lambda) {
    // All shells collapse to R=0 at time tc -- minimized integrand
    // in my written notes
    double a;
    a = pow(3./(2.*lambda),1./3.);
    return a*z;
}

inline double E2_init_focus_rprime(double z, double lambda) {
    // Rprime(0,z) corresponding to above function
    return pow((3./(2*lambda)),1./3.);

}

inline mp_type Shellfocus_sing(mp_type z, mp_type p) {
    // p is a shitty placeholder for the root finder
    mp_type lambda=1.0;
    mp_type m;
    m = pow(z,3.)/2.;
    return lambda*m*m + 3.*pow((mp_type)energy_e2((double)z),3.);
}



int main(int argc, char** argv) {

    // input file streams
    std::ifstream r_init("Example1/R_init.dat");

    // output file streams
    std::ofstream test_num_sol,test2, rdot_f,tsol_f;
    std::ofstream eta_bang_f,eta_crunch_f;
    std::ofstream alan_roots;
    alan_roots.open("alanroots.dat");
    
    std::vector<std::vector<double> > file_input;
    file_input = read_file(r_init);

    test_num_sol.open("test.dat");
    test2.open("test2.dat");
    rdot_f.open("rdot.dat");
    tsol_f.open("tsol.dat");


    /***********************************************************************************
     * DETERMINING INITIAL CONDITIONS FOR INTEGRATION
     * 
     * Solving the exact system at t=0 for example 1
    ************************************************************************************/
    std::vector<std::vector<double> > init_conds_zeta, initial_data, domain;
    std::vector<double> r_init_conds, data_slice;


    // Best way is the compute the zero set directly inside the program. This may take a while.
    // init_conds = compute_zero_set_func(eta_eq);

    // for (size_t i = 0; i < 10; i++)
    // {
    //     r_init_conds.push_back(r_exact(file_input[0][i],file_input[1][i]));
    //     // std::cout << file_input[0][i] << " " << file_input[1][i] << " " << r_init_conds[i] << std::endl;
    //     data_slice.push_back(file_input[0][i]);
    //     data_slice.push_back(r_init_conds[i]);
    //     initial_data.push_back(data_slice);
    //     data_slice.clear();
    // }

    std::ofstream z_out;
    double z_init, z_end, eta_init, eta_end, eta_current, x_grid,y_grid;
    size_t z_n, eta_n;
    std::vector<double> z_grid,eta_grid;
    std::vector<mp_type> z_roots, eta_roots;
    std::vector<std::vector<mp_type> > coarse_domain;
    // grid definitions:
    z_init = 0. +  1e-4;
    z_end = 1.;
    eta_init = 0.7;
    eta_end = 1.2;
    z_n = 60;
    eta_n = 80;

    // Choosing x,y position in the spacetime for the quasispherical case.
    x_grid = 0;
    y_grid = 0;
    /*

        INITIAL CONDITIONS ARE A CRITICAL POINT FOR NON ZERO LAMBDA PROBLEM!
        I.e. SOLUTION WILL NOT GO ANYWHERE

        Need better initial conditions, perhaps theres a better way to pick.

    */

    z_roots = create_grid((mp_type) z_init, (mp_type) z_end, z_n);
    eta_roots = create_grid((mp_type) eta_init, (mp_type) eta_end, eta_n);
    // Create mesh grids:
    z_grid = create_grid(z_init,z_end,z_n);
    // for (int i = 0; i < 200; ++i)
    // {
    //     z_grid.push_back(0.63 + ((double)i/200.)*(0.005));
    // }
    // for (int i = 0; i < 20; ++i)
    // {
    //     z_grid.push_back(0.635 + ((double)i/20.)*(0.165));
    // }

    eta_grid = create_grid(eta_init, eta_end, eta_n);

    coarse_domain.push_back(z_roots);
    coarse_domain.push_back(eta_roots);
    z_out.open("z_out.dat");
    std::ofstream r1_f;
    r1_f.open("debnath.dat");
    for (int i = 0; i < z_grid.size(); ++i)
    {
        z_out << z_grid[i] << std::endl;
    }

    // compute_and_save_zero_set(Shellfocus_sing,coarse_domain,"para.dat");
    // compute_and_save_zero_set(alan_theta,coarse_domain,"para.dat");
    // compute_and_save_zero_set(plane_para_ellipse,time_example1,"para2.dat");
    // compute_and_save_zero_set(alan_theta,coarse_domain,"alan_f.dat");
    // compute_and_save_zero_set(R_example1,time_example1,"alantheta.dat");

    // compute_and_save_zero_set(gaussian,coarse_domain,"alan_f.dat");


    /**
     * ODE INTEGRATION
     * 
     * Solver works on example 1 and example 2, initial conditions must be modified
     * 
    */
    // Integration method
    boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // Solution containers
    std::vector<double> r_curve(2);
    std::vector<double> t_sol;
    std::vector<std::vector<double> > R_sol, full_sol_transformed,Rprime_sol, rdot_sol, fst_refine;
    std::vector<std::vector<std::vector<double> > > full_solution, time_solution, fs_refine;
    std::vector<double> rdot_slice;
    // Initial Conditions and parameter values
    double t_start = 0;
    double t_end = 20;
    double dt = 0.001;
    double lambda =0.1;// pow(2./(3.*pow(1.0,3)),2);
    t_sol = create_grid(t_start,t_end,(int)(t_end-t_start)/dt);

    for (int i = 0; i < t_sol.size(); ++i)
    {
        tsol_f << t_sol[i] << std::endl;
    }

    // Testing the horizon locations with the table in Debnath et al 2006
    for (int i = 0; i < z_grid.size(); ++i)
    {
        r1_f << R2_debnath(z_grid[i],lambda) << " " << z_grid[i] << " " << 2*mass_e2(z_grid[i]) << " " << energy_e2(z_grid[i]) << std::endl;
        r1_f << R1_debnath(z_grid[i],lambda) << " " << z_grid[i] << " " << 2*mass_e2(z_grid[i]) << " " << energy_e2(z_grid[i]) << std::endl;
    }

    std::ofstream t_ah;
    t_ah.open("t_ah.dat");
    int num_R = 1000;
    double R0 = 0;
    double RN = R2_debnath(z_grid[10],lambda);
    double Rtemp;
    std::vector<double> R_domain;
  //   for (int i = 0; i < num_R; ++i)
  //   {
		// Rtemp = R0 + i*(RN-R0)/num_R;
  //   	t_ah << integrand(Rtemp,z_grid[10],lambda) << std::endl;
  //   }
    double z_value;
    std::cout << "Calculating solution curve for i= " << std::endl;
    // Looping through all initial conditions to get a series of solution curves in the z space.
    for (size_t i = 0; i < file_input.size(); i++)
    {
        // Initial conditions for example 1:
        r_curve[0] = file_input[i][1];
        z_value = file_input[i][0];
        std::cout << i << "/" << z_grid.size() << "\r" << std::flush;

        // Initial conditions for example 2: influenced by exact soln
        // r_curve[0] = example2_Rmax(z_grid[i]);
        // r_curve[1] = example2_Rprime_init(z_grid[i]);

        // Debnath and Nolan Comoving frame choice R(0,r)=r, R' = 1
        // r_curve[0] = pow(z_grid[i],2)/2.;
        // r_curve[1] = z_grid[i];

        // Initial conditions to focus R=0 at a single time instant for all shells ? Maybe
        // r_curve[0] = E2_init_focus_r(z_grid[i],lambda);
        // r_curve[1] = E2_init_focus_rprime(z_grid[i],lambda);

        // Compute the collapse time for each shell and use it to limit the integrator
        // std::cout << i << ", t_c = " << trapz(r_curve[0],z_grid[i],lambda) << ", z_i = " << z_grid[i] << std::endl;
        // t_end = trapz(r_curve[0],z_grid[i],lambda);



        // r_curve[0] = zero_force_init(z_grid[i],lambda);
        // r_curve[0] = 10*z_grid[i];

        // ODE to solve for each initial condition. PBH Formation!
        // ode_e2 primordial_bh(z_grid[i],lambda,false,x_grid,y_grid);
        // boost::numeric::odeint::integrate_const(stepper,primordial_bh, r_curve, t_start,t_end,dt,push_back_state_and_time(R_sol,t_sol));
        // full_solution.push_back(R_sol);
        
        // ODE Example 1 GBH Formation.
        ode_e1 galactic_bh(z_value,lambda,false);
        boost::numeric::odeint::integrate_const(stepper,galactic_bh, r_curve, t_start,t_end,dt,push_back_state_and_time(R_sol,t_sol));
        full_solution.push_back(R_sol);


        R_sol.clear();


    }
    std::cout << "Completed ODE solutions" << std::endl;
    std::cout << "Saving solution curves to file" << std::endl;
    // Removing the NaN values from the solution curve vectors.    
    full_sol_transformed = removeNAN(full_solution);
    Rprime_sol = removeNAN(full_solution,1);
    matrix_to_file3(full_sol_transformed,"test2.dat");
    matrix_to_file3(Rprime_sol,"Rprime.dat");

    /*
        OUTPUT OF THE INVARIANTS AND THE APPARENT HORIZON EQUATION

    */

    std::cout << "Calculating Horizon and Shell Crossing Detectors" << std::endl;
    std::cout << "Check the functions used across examples!!" << std::endl;

    std::vector<std::vector<double> > mu_sol, rho_sol,app_sol;
    std::vector<double> mu_slice, rho_slice, app_slice;
    double rdot_i, e_i, r_i,mu_i,app_i,z_i;
    // Output some extra things like rho and mu
    std::cout << full_sol_transformed.size() << std::endl;
    for (int i = 0; i < full_sol_transformed.size(); ++i)
    {

        for (int j = 0; j < full_sol_transformed[i].size(); ++j)
        {
            z_i = file_input[i][0];
            r_i = full_sol_transformed[i][j];
            rdot_i = Rdot_e1(r_i,z_i,lambda);
            

            // Change energy_e2 to energy_e1 for proper analysis.
            e_i = energy_e11(z_i);
            mu_i = mu_example2(rdot_i,r_i,e_i);
            // mu_i = apparent_horizon(r_i,file_input[i][0],lambda);
            // rdot_slice.push_back(Rdot(full_sol_transformed[i][j],file_input[i][0],lambda));
            mu_slice.push_back(mu_i);

            app_i = apparent_horizon_e1(r_i,z_i,lambda);
            app_slice.push_back(app_i);
        }
        mu_sol.push_back(mu_slice);
        mu_slice.clear();
        app_sol.push_back(app_slice);
        app_slice.clear();
    }
    zeros_output(full_sol_transformed, mu_sol,z_grid,t_sol,"mu_zeros.dat",lambda);
    zeros_output(full_sol_transformed, app_sol,z_grid,t_sol,"apparent_zeros.dat",lambda);
    matrix_to_file3(mu_sol,"rdot.dat");

     /*
       Finding the zeros of the Rprime solution curve, this is C_4 !!

    */
    // Rescale Y' solution to (R' - RH'/H) = 0 for proper quasispherical shell crossing
    // std::vector<double> Rprime_slice;
    // std::vector<std::vector<double> > Rprime_rescaled;
    // double z_temp, H_temp,Rprime_temp;

    // for (int i = 0; i < Rprime_sol.size(); ++i)
    // {
    //     z_temp = z_grid[i];
    //     H_temp = H(z_temp,x_grid,y_grid);
    //     for (int j = 0; j < full_sol_transformed[i].size()-1; ++j)
    //     {
    //         // std::cout << "i,j = " << i << " " << j << std::endl; 
    //         Rprime_temp = Rprime_sol[i][j] - full_sol_transformed[i][j]*Hprime(z_temp,x_grid,y_grid)/H(z_temp,x_grid,y_grid);
    //         Rprime_slice.push_back(Rprime_temp);

    //     }
    //     Rprime_rescaled.push_back(Rprime_slice);
    //     Rprime_slice.clear();
    // }

    // matrix_to_file3(Rprime_rescaled,"Yprime.dat");
    // zeros_output(full_sol_transformed, Rprime_sol,z_grid,t_sol,"Rprime_zeros.dat",lambda);
    // zeros_output(full_sol_transformed, Rprime_rescaled, z_grid, t_sol, "Rprime_rescaled_zeros.dat",lambda);





    /*
		Refine the bracketed roots given above using the Horizon detectors(mu=0, AH=0) 
		by evolving the differential equation again, starting at R_leftbracket(z,t_leftbracket)
		and evolving up until the point in time t_rightbracket, and say a large number of steps 
		in between to determine a finer bracket for mu and AH.
    */


 //    double t_start_refine, t_end_refine, dt_refine;
 //    int N_refine;
 //    N_refine = 100;
 //    for (int i = 0; i < mu_sol.size(); ++i)
 //    {
 //    	// Set up the ODE again
 //    	ode_e2 Refiningroots(mu_sol[i][1],lambda,false);

	// 	// Define initial conditions. Starting at R in mu_sol
 //    	r_curve[0] = mu_sol[i][0];
 //    	t_start_refine = mu_sol[i][2];
 //    	t_end_refine = mu_sol[i][3];
 //    	dt_refine = (t_end_refine - t_start_refine)/N_refine;

 //    	boost::numeric::odeint::integrate_const(stepper,Refiningroots, r_curve, t_start_refine,t_end_refine,dt_refine,push_back_state_and_time(R_sol,t_sol));

 //    	fs_refine.push_back(R_sol);
 //    	R_sol.clear();
 //    }

 //    fst_refine = removeNAN(fs_refine);
 //    matrix_to_file3(fst_refine, "mu_refined_sol.dat");


 //    std::vector<std::vector<double> > mu_sol_r, rho_sol_r,app_sol_r;
 //    std::vector<double> mu_slice_r, rho_slice_r, app_slice_r;


	// for (int i = 0; i < fst_refine.size(); ++i)
 //    {

 //        for (int j = 0; j < fst_refine[i].size(); ++j)
 //        {
 //            rdot_i = Rdot(fst_refine[i][j],mu_sol[i][1],lambda);
 //            r_i = fst_refine[i][j];
 //            e_i = energy_e2(mu_sol[i][1]);
 //            mu_i = mu_example2(rdot_i,r_i,e_i);
 //            // mu_i = apparent_horizon(r_i,mu_sol[i][1],lambda);
 //            // rdot_slice.push_back(Rdot(full_sol_transformed[i][j],mu_sol[i][1],lambda));
 //            mu_slice_r.push_back(mu_i);

 //            app_i = apparent_horizon(r_i,mu_sol[i][1],lambda);
 //            app_slice_r.push_back(app_i);
 //        }
 //        mu_sol_r.push_back(mu_slice_r);
 //        mu_slice_r.clear();
 //        app_sol_r.push_back(app_slice_r);
 //        app_slice_r.clear();
 //    }

 //    zeros_output(fst_refine, mu_sol_r,z_grid,t_sol,"mu_zeros_R.dat");
 //    zeros_output(fst_refine, app_sol_r,z_grid,t_sol,"apparent_zeros_R.dat");







} // End Main


void compute_and_save_zero_set(mp_type (*f)(mp_type,mp_type), std::vector<std::vector<mp_type> > domain, std::string filename) {
    // Input function: f = parameterization of R
    // Input domain: 2xN vector of meshgrid on which to compute the zero set.                 
    // Mesh in (z,eta) space
    std::vector<std::vector<mp_type> > zero_set;
    std::ofstream outfile;
    outfile.open(filename);

    // Measuring compute time for the root finder.
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "Bracketing roots \n";
    // bracket the roots of the surface
    zero_set = bracket_secant_method(f,domain[0],domain[1],6);


    auto current_time = std::chrono::high_resolution_clock::now();
    std::cout << "# Writing results to file" << std::endl;
    std::cout << "# X domain: [" << domain[0][0] << "," << domain[0][domain[0].size()-1] << "]\n";
    std::cout << "# Y domain: [" << domain[1][0] << "," << domain[1][domain[1].size()-1] << "]\n";
    std::cout << "# Grid size (X x Y): " << domain[0].size() << " x " << domain[1].size() << std::endl;
    std::cout << "# Roots found: " << zero_set[0].size() << std::endl;
    std::cout << "# Compute time: " << std::chrono::duration_cast<std::chrono::milliseconds>((current_time - start_time)).count() << " ms" << std::endl;
    
    outfile << "# Writing results to file" << std::endl;
    outfile << "# X domain: [" << domain[0][0] << "," << domain[0][domain[0].size()-1] << "]\n";
    outfile << "# Y domain: [" << domain[1][0] << "," << domain[1][domain[1].size()-1] << "]\n";
    outfile << "# Grid size (X x Y): " << domain[0].size() << " x " << domain[1].size() << std::endl;
    outfile << "# Roots found: " << zero_set[0].size() << std::endl;
    outfile << "# Compute time: " << std::chrono::duration_cast<std::chrono::milliseconds>((current_time - start_time)).count() << " ms" << std::endl;
    outfile << "# Data Column layout: {r_0,z_0,r_bl,z_bl,r_br,z_br}\n";
    // output resulting set to file
    matrix_to_file(zero_set,outfile);

}

std::vector<std::vector<double> > read_file(std::ifstream& input) {
    // Important note: This reads across rows, then columns, stores nxn array into nx1.
    // I need to transform the input data back to nx2
    std::istream_iterator<double> start(input), end;
    std::vector<double> numbers(start,end);
    std::vector<double> row;
    std::vector<std::vector<double> > result;
    std::cout << "Stored " << numbers.size() << " values from file." << std::endl;

    // Captures the left column
    for (size_t i = 0; i < numbers.size(); i+=2) {
        row.push_back(numbers[i]);
        row.push_back(numbers[i+1]);
        result.push_back(row);
        row.clear();
    }

    return result;
}


