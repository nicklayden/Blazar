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


inline mp_type alan_f(mp_type r, mp_type k) {
    // Alans specific prob. Just wanted a graph showing the function
    // crossing f(r,k)=0 for r near pi/4 and k = (1,2)
    mp_type a,b,c,d,e,f,g,h;
    a = sinh(k*(M_PI - r) + sinh(k*r));
    b = cosh(k*(M_PI - r) - cosh(k*r));
    c = sinh(k*(M_PI/2. + r) + sinh(k*(M_PI/2. - r)));
    d = cosh(k*(M_PI/2. + r) - cosh(k*(M_PI/2. - r)));

    e = -(1./sin(r));
    f = -cos(r)/(sin(r)*sin(r));
    g = 1./cos(r);
    h = tan(r);

    return e*(  f*a + k*b) + g*( h*c + k*d )  ;
}


inline mp_type alan_f_A(mp_type r, mp_type k) {
    // Alans function
    // THIS IS A(r) 
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


int main(int argc, char** argv) {

    // input file streams
    std::ifstream r_init("R_init.dat");

    // output file streams
    std::ofstream test_num_sol,test2;
    std::ofstream eta_bang_f,eta_crunch_f;
    std::ofstream alan_roots;
    alan_roots.open("alanroots.dat");
    
    std::vector<std::vector<double> > file_input;
    file_input = read_file(r_init);

    test_num_sol.open("test.dat");
    test2.open("test2.dat");
    /**
     * DETERMINING INITIAL CONDITIONS FOR INTEGRATION
     * 
     * Solving the exact system at t=0 for example 1
    */
    std::vector<std::vector<mp_type> > init_conds_zeta, initial_data, domain;
    std::vector<mp_type> r_init_conds, data_slice;



    // init_conds = compute_zero_set_func(eta_eq);

    for (size_t i = 0; i < 10; i++)
    {
        r_init_conds.push_back(r_exact(file_input[0][i],file_input[1][i]));
        // std::cout << file_input[0][i] << " " << file_input[1][i] << " " << r_init_conds[i] << std::endl;
        data_slice.push_back(file_input[0][i]);
        data_slice.push_back(r_init_conds[i]);
        initial_data.push_back(data_slice);
        data_slice.clear();
    }


    mp_type z_init, z_end, eta_init, eta_end, eta_current;
    size_t z_n, eta_n;
    std::vector<mp_type> z_grid,eta_grid;
    std::vector<std::vector<mp_type> > coarse_domain;
    // grid definitions:
    z_init = -5;
    z_end = 5;
    eta_init = -5;
    eta_end = 5;
    z_n = 100;
    eta_n = 100;

    // Create mesh grids:
    z_grid = create_grid(z_init,z_end,z_n);
    eta_grid = create_grid(eta_init, eta_end, eta_n);

    coarse_domain.push_back(z_grid);
    coarse_domain.push_back(eta_grid);

    // compute_and_save_zero_set(test_paraboloid,time_example1,"para.dat");
    // compute_and_save_zero_set(plane_para_ellipse,time_example1,"para2.dat");
    // compute_and_save_zero_set(alan_ar_test,time_example1,"alan_f.dat");
    // compute_and_save_zero_set(R_example1,time_example1,"alantheta.dat");

    compute_and_save_zero_set(gaussian,coarse_domain,"alan_f.dat");


    // /**
    //  * ODE INTEGRATION
    //  * 
    //  * First for Example 1 in Szekeres Paper
    //  * 
    // */
    // // Integration method
    // boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // // Solution containers
    // std::vector<double> r_curve(1);
    // std::vector<double> t_sol;
    // std::vector<std::vector<double> > R_sol, full_sol_transformed;
    // std::vector<std::vector<std::vector<double> > > full_solution;
    // // Initial Conditions
    // double t_start = 0;
    // double t_end = 18;
    // double dt = 0.001;
    // // int kk = 483;
    // // Integrate through all z, each solution is unique for each z.
    // for (size_t i = 0; i < file_input.size(); i++)
    // {
    //     // Set initial conditions for R and z.
    //     r_curve[0] = file_input[i][1];
    //     // std::cout << r_curve[0] << std::endl;
    //     ode_e1 solution_curve(file_input[i][0],false);
    //     boost::numeric::odeint::integrate_const(stepper,solution_curve, r_curve, t_start,t_end,dt,push_back_state_and_time(R_sol,t_sol));
    //     // The above line gives a solution curve for a single z value. We need to iterate
    //     // through the grid of z values to complete the curve.
    //     // NOTE: Need to scroll through solution curves and remove all NAN values.
    //     full_solution.push_back(R_sol);
    //     R_sol.clear();

    // }
    
    // std::cout << full_solution[0][1][0] << isnan(full_solution[0][1][0]) <<  std::endl;
    // if (!isnan(full_solution[0][0][0]))
    // {
    //     std::cout << full_solution[0][0][0] << std::endl;   
    // }

    // // Removing the NaN values from the solution curve vectors.    
    // full_sol_transformed = removeNAN(full_solution);

    // // lbracket,rbracket is R(z) for either side of the root
    // // zbracket, tbracket are the corresponding z and t values for the lbracket
    // std::vector<double> lbracket,rbracket, zbracket, tbracket;
    // double lvalue,rvalue;
    // // Bracketing the solution sections to determine the solution for R=2M
    // for (size_t i = 0; i < full_sol_transformed.size(); i++) {
    //     if (full_sol_transformed[i].size() > 1) {
    //         for (size_t j = 1; j < full_sol_transformed[i].size(); j++) {
    //             lvalue = full_sol_transformed[i][j-1] - 2.* file_input[i][0];
    //             rvalue = full_sol_transformed[i][j] - 2.* file_input[i][0];
    //             if (lvalue*rvalue < 0.) {
    //                 lbracket.push_back(full_sol_transformed[i][j-1]);
    //                 rbracket.push_back(full_sol_transformed[i][j]);
    //                 zbracket.push_back(file_input[i][0]);
    //                 tbracket.push_back(t_sol[j]);
    //             }
    //         }
    //     }
    // }
    
    // // std::cout << lbracket.size() << " " << rbracket.size() << " " << zbracket.size() << " " << tbracket.size() << std::endl;


    // // for (size_t i = 0; i < lbracket.size(); i++)
    // // {
    // //     std::cout << lbracket[i] << " " << rbracket[i] << " " << zbracket[i] << " " << tbracket[i] << " " << std::endl;
    // // }
    

    // matrix_to_file(full_sol_transformed,test2);

    // // for ( size_t i = 0; i < 30; i++) { 
    // //     for ( size_t j = 0; j < full_sol_transformed[i].size(); j++) {
    // //         std::cout << full_sol_transformed[i][j] << " ";
    // //     }
    // //     std::cout << std::endl;
    // // }

    // std::cout << full_sol_transformed.size() << " " << full_sol_transformed[0].size() << " " << full_sol_transformed[100].size() << std::endl;


    // std::cout << full_solution.size() << " " << full_solution[0].size() << std::endl;

    // // matrix_to_file2(full_solution,test_num_sol);

    // // for (size_t i = 0; i < R_sol.size(); i++)
    // // {
    // //     std::cout << R_sol[i][0] << " " << R_sol[i][0] - 2.*file_input[kk][0] << " " << t_sol[i] <<  std::endl;
    // //     // for (size_t j = 0; j < R_sol[0].size(); j++)
    // //     // {
    // //     //     std::cout << R_sol[i][0] << " ";
    // //     // }
    // //     // std::cout << std::endl;
        
    // // }
    







    // for (size_t i = 0; i < file_input.size(); i++)
    // {
    //     for (size_t j = 0; j < file_input[0].size(); j++)
    //     {
    //         std::cout << file_input[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    

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
    zero_set = bracket_secant_method(f,domain[0],domain[1],4);


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






