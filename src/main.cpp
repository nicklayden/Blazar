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

inline double Rdot(double r, double z, double lambda) {
    double a;
    a = 2*energy_e2(z) + 2*mass_e2(z)/r + lambda*r*r/3.;
    return -sqrt(a);
}

double mu_example2(double rdot, double r, double energy) {
    // Extended Cartan invariant that detects the horizon
    double a,b,c;
    a = sqrt(1 + 2*energy);
    b = rdot + a;
    c = b/r;
    return -c/sqrt(2.0);
}

double apparent_horizon(double R, double z, double lambda) {
    // Debnath Nath Chakraborty 2006 paper equaiton 25
    // Requires mass for EXAMPLE 2
    return lambda*pow(R,3) + 6*mass_e2(z) - 3*R;


}

void zeros_output(std::vector<std::vector<double> > rsol, std::vector<std::vector<double> > function,std::vector<double> z, std::vector<double> t, std::string file) {
    // Compute simple bracketing set , return one half as a pseudo zero set
    std::vector<std::vector<double> > zeros;
    std::vector<double> slice;
    double lhs,rhs;

    for (int i = 0; i < function.size(); ++i)
    {
        slice.clear();
        for (int j = 0; j < function[i].size()-1; ++j)
        {
            lhs = function[i][j];
            rhs = function[i][j+1];
            if (lhs * rhs <= 0.0) {
                slice.push_back(rsol[i][j]);
                slice.push_back(z[i]);
                slice.push_back(t[j]);
                zeros.push_back(slice);
                slice.clear();
            }
            

        }
        // zeros.push_back(slice);
        // slice.clear();
    }
    matrix_to_file3(zeros,file);



}


double R1_debnath(double z, double lambda) {
    double a,b,c;
    a = 2./sqrt(lambda);

    c = acos(-(3./2.)*sqrt(lambda*2*mass_e2(z)));
    b = cos(c/3.);
    return a*b;
}

double R2_debnath(double z, double lambda) {
    double a,b,c,theta;
    a = 1./sqrt(lambda);
    theta = acos(-(3./2.)*sqrt(lambda*2*mass_e2(z)));
    b = -cos(theta/3.) + sqrt(3.) * sin(theta/3.);
    return a*b;
}





int main(int argc, char** argv) {

    // input file streams
    std::ifstream r_init("R_init.dat");

    // output file streams
    std::ofstream test_num_sol,test2, rdot_f,tsol_f;
    std::ofstream eta_bang_f,eta_crunch_f;
    std::ofstream alan_roots;
    alan_roots.open("alanroots.dat");
    
    // std::vector<std::vector<double> > file_input;
    // file_input = read_file(r_init);

    test_num_sol.open("test.dat");
    test2.open("test2.dat");
    rdot_f.open("rdot.dat");
    tsol_f.open("tsol.dat");


    /***********************************************************************************
     * DETERMINING INITIAL CONDITIONS FOR INTEGRATION
     * 
     * Solving the exact system at t=0 for example 1
    ************************************************************************************/
    std::vector<std::vector<mp_type> > init_conds_zeta, initial_data, domain;
    std::vector<mp_type> r_init_conds, data_slice;



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
    double z_init, z_end, eta_init, eta_end, eta_current;
    size_t z_n, eta_n;
    std::vector<double> z_grid,eta_grid;
    std::vector<std::vector<double> > coarse_domain;
    // grid definitions:
    z_init = 0. +  1e-4;
    z_end = 0.88;
    eta_init = 3.;
    eta_end = 5.;
    z_n = 200;
    eta_n = 100;


    /*

        INITIAL CONDITIONS ARE A CRITICAL POINT FOR NON ZERO LAMBDA PROBLEM!
        I.e. SOLUTION WILL NOT GO ANYWHERE

        Need better initial conditions, perhaps theres a better way to pick.

    */



    // Create mesh grids:
    z_grid = create_grid(z_init,z_end,z_n);
    eta_grid = create_grid(eta_init, eta_end, eta_n);

    coarse_domain.push_back(z_grid);
    coarse_domain.push_back(eta_grid);
    z_out.open("z_out.dat");
    std::ofstream r1_f;
    r1_f.open("debnath.dat");
    for (int i = 0; i < z_grid.size(); ++i)
    {
        z_out << z_grid[i] << std::endl;
    }


    // compute_and_save_zero_set(alan_theta,coarse_domain,"para.dat");
    // compute_and_save_zero_set(plane_para_ellipse,time_example1,"para2.dat");
    // compute_and_save_zero_set(alan_ar_test,time_example1,"alan_f.dat");
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
    std::vector<double> r_curve(1);
    std::vector<double> t_sol;
    std::vector<std::vector<double> > R_sol, full_sol_transformed, rdot_sol;
    std::vector<std::vector<std::vector<double> > > full_solution, time_solution;
    std::vector<double> rdot_slice;
    // Initial Conditions and parameter values
    double t_start = 0;
    double t_end = 20;
    double dt = 0.0001;
    double lambda = 1.0;
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


    // Looping through all initial conditions to get a series of solution curves in the z space.
    for (size_t i = 0; i < z_grid.size(); i++)
    {
        // Initial conditions for example 1:
        // r_curve[0] = file_input[i][1];

        // Initial conditions for example 2:
        r_curve[0] = example2_Rmax(z_grid[i]);
        // r_curve[0] = 10*z_grid[i];

        // ODE to solve for each initial condition.
        ode_e2 primordial_bh(z_grid[i],lambda,false);
        boost::numeric::odeint::integrate_const(stepper,primordial_bh, r_curve, t_start,t_end,dt,push_back_state_and_time(R_sol,t_sol));
        full_solution.push_back(R_sol);
        // std::cout << R_sol[0].size() << std::endl;
        R_sol.clear();

    }

    // Removing the NaN values from the solution curve vectors.    
    full_sol_transformed = removeNAN(full_solution);
    matrix_to_file3(full_sol_transformed,"test2.dat");
    



    /*
        OUTPUT OF THE INVARIANTS AND THE APPARENT HORIZON EQUATION

    */
    std::vector<std::vector<double> > mu_sol, rho_sol,app_sol;
    std::vector<double> mu_slice, rho_slice, app_slice;
    double rdot_i, e_i, r_i,mu_i,app_i;
    // Output some extra things like rho and mu
    for (int i = 0; i < full_sol_transformed.size(); ++i)
    {

        for (int j = 0; j < full_sol_transformed[i].size(); ++j)
        {
            rdot_i = Rdot(full_sol_transformed[i][j],z_grid[i],lambda);
            r_i = full_sol_transformed[i][j];
            e_i = energy_e2(z_grid[i]);
            mu_i = mu_example2(rdot_i,r_i,e_i);
            // mu_i = apparent_horizon(r_i,z_grid[i],lambda);
            // rdot_slice.push_back(Rdot(full_sol_transformed[i][j],z_grid[i],lambda));
            mu_slice.push_back(mu_i);

            app_i = apparent_horizon(r_i,z_grid[i],lambda);
            app_slice.push_back(app_i);
        }
        mu_sol.push_back(mu_slice);
        mu_slice.clear();
        app_sol.push_back(app_slice);
        app_slice.clear();
    }
    zeros_output(full_sol_transformed, mu_sol,z_grid,t_sol,"mu_zeros.dat");
    zeros_output(full_sol_transformed, app_sol,z_grid,t_sol,"apparent_zeros.dat");
    matrix_to_file3(mu_sol,"rdot.dat");

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






