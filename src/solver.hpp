#pragma once

/*
	ODE Solver custom integration function utilizing Boost::odeint integrators


*/


// Boost libraries
#include <boost/multiprecision/mpfr.hpp>
// Note: in mpfr, expression templates are disabled. This reduces compiles but increases runtime.
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint.hpp>

// Standard Libraries
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>

#include "odes.hpp"
#include "methods.hpp"


void Solve() {

	std::vector<double> z_grid;


    // Integration method
    boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // System
	ode_e2 primordial_bh(z_grid[i],lambda,false);




}







