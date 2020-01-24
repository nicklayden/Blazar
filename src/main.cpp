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

// Type definitions
// typedef boost::multiprecision::cpp_dec_float_50 mp_type;
namespace mp = boost::multiprecision;
typedef mp::mpfr_float_100 mp_type; // faster compile time


// Preprocessor defines
#define PI M_PIl

void compute_and_save_zero_set(mp_type (*f)(mp_type,mp_type), mp_type (*g)(mp_type, mp_type), std::ofstream& outfile);
std::vector<std::vector<mp_type> > compute_zero_set_func(mp_type (*f)(mp_type,mp_type));
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

    return e*(  f*a + k*b + g*( h*c + k*d )  );
}



int main(int argc, char** argv) {
    // Set multiprecision digits
    // std::setprecision(std::numeric_limits<mp_type>::max_digits10);
    // mpf_set_default_prec(50);

    // input file streams
    std::ifstream r_init("R_init.dat");

    // output file streams
    std::ofstream rhop_f,rhom_f,psim_f,psip_f,app_et_f,app_ct_f,e1_R,rho_ex1_f,mu_ex1_f,eta_set_ex1,test_num_sol,test2;
    std::ofstream eta_bang_f,eta_crunch_f;
    std::ofstream alan_roots;
    alan_roots.open("alanroots.dat");

    // eta_bang_f.open("eta_bang.dat");
    // eta_crunch_f.open("eta_crunch.dat");

    // rhop_f.open("rhop.dat");
    // rhom_f.open("rhom.dat");
    // psim_f.open("psim.dat");
    // psip_f.open("psip.dat");
    // app_et_f.open("app_et.dat");
    // app_ct_f.open("app_ct.dat");
    // e1_R.open("e1_R.dat");
    // rho_ex1_f.open("rho_ex1.dat");
    // mu_ex1_f.open("mu_ex1.dat");
    // eta_set_ex1.open("eta_set_ex1.dat");
    std::vector<std::vector<double> > file_input;
    file_input = read_file(r_init);
    std::cout << "Sublime test\n";
    test_num_sol.open("test.dat");
    test2.open("test2.dat");
    /**
     * DETERMINING INITIAL CONDITIONS FOR INTEGRATION
     * 
     * Solving the exact system at t=0 for example 1
    */
    std::vector<std::vector<mp_type> > init_conds_zeta, initial_data;
    std::vector<mp_type> r_init_conds, data_slice;

    // init_conds = compute_zero_set_func(eta_eq);

    for (size_t i = 0; i < 10; i++)
    {
        r_init_conds.push_back(r_exact(file_input[0][i],file_input[1][i]));
        std::cout << file_input[0][i] << " " << file_input[1][i] << " " << r_init_conds[i] << std::endl;
        data_slice.push_back(file_input[0][i]);
        data_slice.push_back(r_init_conds[i]);
        initial_data.push_back(data_slice);
        data_slice.clear();
    }
    
    compute_and_save_zero_set(alan_f,time_example1,alan_roots);
    // compute_and_save_zero_set(eta_bang,time_example1,eta_bang_f);
    // compute_and_save_zero_set(eta_crunch,time_example1,eta_crunch_f);    
    // compute_and_save_zero_set(eta_eq,time_example1,eta_set_ex1);


    /**
     * ODE INTEGRATION
     * 
     * First for Example 1 in Szekeres Paper
     * 
    */
    // Integration method
    boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // Solution containers
    std::vector<double> r_curve(1);
    std::vector<double> t_sol;
    std::vector<std::vector<double> > R_sol, full_sol_transformed;
    std::vector<std::vector<std::vector<double> > > full_solution;
    // Initial Conditions
    double t_start = 0;
    double t_end = 18;
    double dt = 0.001;
    // int kk = 483;
    // Integrate through all z, each solution is unique for each z.
    for (size_t i = 0; i < file_input.size(); i++)
    {
        // Set initial conditions for R and z.
        r_curve[0] = file_input[i][1];
        // std::cout << r_curve[0] << std::endl;
        ode_e1 solution_curve(file_input[i][0],false);
        boost::numeric::odeint::integrate_const(stepper,solution_curve, r_curve, t_start,t_end,dt,push_back_state_and_time(R_sol,t_sol));
        // The above line gives a solution curve for a single z value. We need to iterate
        // through the grid of z values to complete the curve.
        // NOTE: Need to scroll through solution curves and remove all NAN values.
        full_solution.push_back(R_sol);
        R_sol.clear();

    }
    
    std::cout << full_solution[0][1][0] << isnan(full_solution[0][1][0]) <<  std::endl;
    if (!isnan(full_solution[0][0][0]))
    {
        std::cout << full_solution[0][0][0] << std::endl;   
    }

    // Removing the NaN values from the solution curve vectors.    
    full_sol_transformed = removeNAN(full_solution);

    // lbracket,rbracket is R(z) for either side of the root
    // zbracket, tbracket are the corresponding z and t values for the lbracket
    std::vector<double> lbracket,rbracket, zbracket, tbracket;
    double lvalue,rvalue;
    // Bracketing the solution sections to determine the solution for R=2M
    for (size_t i = 0; i < full_sol_transformed.size(); i++) {
        if (full_sol_transformed[i].size() > 1) {
            for (size_t j = 1; j < full_sol_transformed[i].size(); j++) {
                lvalue = full_sol_transformed[i][j-1] - 2.* file_input[i][0];
                rvalue = full_sol_transformed[i][j] - 2.* file_input[i][0];
                if (lvalue*rvalue < 0.) {
                    lbracket.push_back(full_sol_transformed[i][j-1]);
                    rbracket.push_back(full_sol_transformed[i][j]);
                    zbracket.push_back(file_input[i][0]);
                    tbracket.push_back(t_sol[j]);
                }
            }
        }
    }
    
    // std::cout << lbracket.size() << " " << rbracket.size() << " " << zbracket.size() << " " << tbracket.size() << std::endl;


    // for (size_t i = 0; i < lbracket.size(); i++)
    // {
    //     std::cout << lbracket[i] << " " << rbracket[i] << " " << zbracket[i] << " " << tbracket[i] << " " << std::endl;
    // }
    

    matrix_to_file(full_sol_transformed,test2);

    // for ( size_t i = 0; i < 30; i++) { 
    //     for ( size_t j = 0; j < full_sol_transformed[i].size(); j++) {
    //         std::cout << full_sol_transformed[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::cout << full_sol_transformed.size() << " " << full_sol_transformed[0].size() << " " << full_sol_transformed[100].size() << std::endl;


    std::cout << full_solution.size() << " " << full_solution[0].size() << std::endl;

    // matrix_to_file2(full_solution,test_num_sol);

    // for (size_t i = 0; i < R_sol.size(); i++)
    // {
    //     std::cout << R_sol[i][0] << " " << R_sol[i][0] - 2.*file_input[kk][0] << " " << t_sol[i] <<  std::endl;
    //     // for (size_t j = 0; j < R_sol[0].size(); j++)
    //     // {
    //     //     std::cout << R_sol[i][0] << " ";
    //     // }
    //     // std::cout << std::endl;
        
    // }
    







    // for (size_t i = 0; i < file_input.size(); i++)
    // {
    //     for (size_t j = 0; j < file_input[0].size(); j++)
    //     {
    //         std::cout << file_input[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    


    // compute_and_save_zero_set(mu_example1,time_example1,mu_ex1_f);
    // compute_and_save_zero_set(R_example1,time_example1,e1_R);
    // compute_and_save_zero_set(eta_eq,time_example1,eta_set_ex1);
    // compute_and_save_zero_set(rho_ex1,time_example1,rho_ex1_f);
    // compute_and_save_zero_set(psim,ct_function,psim_f);
    // compute_and_save_zero_set(psip,et_function,psip_f);
    // compute_and_save_zero_set(app_horizon,ct_function,app_ct_f);
    // compute_and_save_zero_set(app_horizon,et_function,app_et_f);
    // compute_and_save_zero_set(rhom,ct_function,rhom_f); // note: rhom not defined yet
    // compute_and_save_zero_set(rhop,et_function,rhop_f);

}


void compute_and_save_zero_set(mp_type (*f)(mp_type,mp_type), mp_type (*g)(mp_type, mp_type), std::ofstream& outfile) {
    // Input functions: f = parameterization of R
    //                  g = function from (z,eta)-> (t)
    // Mesh in (z,eta) space
    std::vector<mp_type> z_grid, eta_grid;
    std::vector<std::vector<mp_type> > surface_ze, bracketing_set, zero_set, zero_zt_set;
    auto start_time = std::chrono::high_resolution_clock::now();
   

   
    // Mesh parameters, limits / number of points on grid
    mp_type z_init, z_end, eta_init, eta_end, eta_current;
    size_t z_n, eta_n;

    // grid definitions:
    // z_init = 0.;// + 1e-4;
    // z_end = 5.0;// - 1e-4;
    // eta_init = 0;// + 1e-4;
    // eta_end  = 2*PI - 1e-4;
    z_init = M_PI/4. - 0.1;
    z_end = M_PI/4. + 0.1;
    eta_init = 1.01;
    eta_end = 2.0;
    z_n = 500;
    eta_n = 500;

    // Create mesh grids:
    z_grid = create_grid(z_init,z_end,z_n);
    eta_grid = create_grid(eta_init, eta_end, eta_n);

    std::cout << "Computing surface \n";
    // compute the surfaces
    surface_ze = compute_surface(f,z_grid,eta_grid);

    std::cout << "Bracketing roots \n";
    // bracket the roots of the surface
    bracketing_set = bracket_roots(surface_ze,z_grid,eta_grid);
    // compute the zeros
    zero_set = compute_zero_set(f,bracketing_set,5);

    zero_zt_set.push_back(zero_set[0]);
    zero_zt_set.push_back(zero_set[1]);
    // Transforming from (z,eta) -> (z,t) coordinates (OR NOT - REDO FUNCTION)
    std::cout << "NOTE: NOT TRANSFORMING RESULTING DATA" << std::endl;
    // for (size_t i = 0; i < zero_set[0].size(); i++)
    // {
    //     zero_zt_set[0][i] = zero_zt_set[0][i];
    //     // zero_zt_set[1][i] = g(zero_set[0][i],zero_set[1][i]);
    //     zero_zt_set[1][i] = zero_set[1][i]; // This leaves the coord system alone and DOESNT TRANSFORM!
    //     // std::cout << zero_zt_set[0][i] << " " << zero_zt_set[1][i] << std::endl;
    // }
    


    auto current_time = std::chrono::high_resolution_clock::now();
    std::cout << "Writing results to file" << std::endl;
    std::cout << "Time to compute zero set: " << std::chrono::duration_cast<std::chrono::milliseconds>((current_time - start_time)).count() << " ms" << std::endl;
    // output resulting set to file
    matrix_to_file(zero_zt_set,outfile);

}

std::vector<std::vector<mp_type> > compute_zero_set_func(mp_type (*f)(mp_type,mp_type)) {
    // Input functions: f = parameterization of R
    //                  g = function from (z,eta)-> (t)
    // Mesh in (z,eta) space
    std::vector<mp_type> z_grid, eta_grid;
    std::vector<std::vector<mp_type> > surface_ze, bracketing_set, zero_set, zero_zt_set;

    // Mesh parameters, limits / number of points on grid
    mp_type z_init, z_end, eta_init, eta_end, eta_current;
    size_t z_n, eta_n;

    // grid definitions:
    z_init = 0.;// + 1e-4;
    z_end = 5.0;// - 1e-4;
    eta_init = 0;// + 1e-4;
    eta_end  = 2*PI;// - 1e-4;
    z_n = 500;
    eta_n = 500;

    // Create mesh grids:
    z_grid = create_grid(z_init,z_end,z_n);
    eta_grid = create_grid(eta_init, eta_end, eta_n);


    // compute the surfaces
    surface_ze = compute_surface(f,z_grid,eta_grid);

    // bracket the roots of the surface
    bracketing_set = bracket_roots(surface_ze,z_grid,eta_grid);
    // compute the zeros
    zero_set = compute_zero_set(f,bracketing_set,5);

    zero_zt_set.push_back(zero_set[0]);
    zero_zt_set.push_back(zero_set[1]);
    // Transforming from (z,eta) -> (z,t) coordinates (OR NOT - REDO FUNCTION)
    // std::cout << "NOTE: NOT TRANSFORMING RESULTING DATA" << std::endl;
    // for (size_t i = 0; i < zero_set[0].size(); i++)
    // {
    //     zero_zt_set[0][i] = zero_zt_set[0][i];
    //     // zero_zt_set[1][i] = g(zero_set[0][i],zero_set[1][i]);
    //     zero_zt_set[1][i] = zero_set[1][i]; // This leaves the coord system alone and DOESNT TRANSFORM!
    //     // std::cout << zero_zt_set[0][i] << " " << zero_zt_set[1][i] << std::endl;
    // }
    return zero_zt_set;
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





