#pragma once
/**
 * Functions for the Cartan invariants.
 * Created from the maple output with CodeGeneration + a custom
 * python script
 * 
*/
// Boost libraries
#include <boost/multiprecision/mpfr.hpp>
#include <complex>
// #include <gmpxx.h>
// Type definitions

namespace mp = boost::multiprecision;

// typedef mp::cpp_dec_float_50 mp_type;
typedef mp::mpfr_float_100 mp_type;
typedef std::complex<mp_type> mp_complex;

// typedef boost::multiprecision::cpp_dec_float_100 mp_100;
// typedef mpz_t mp_type;
// R = 2M apparent horizon
mp_type app_horizon(mp_type z, mp_type eta);
mp_type R_example1(mp_type M, mp_type eta);
mp_type energy_e1(mp_type z);
mp_type r_exact(mp_type z, mp_type eta);
mp_type eta_bang(mp_type z, mp_type eta); // equation for eta with the big bang function
mp_type eta_crunch(mp_type z, mp_type eta); // same as above but with the crunch time function

// Cartan Invariants (Contains functions from DIFFERENT problems! Stated ex1 or if not, it belongs to ex 2)
mp_type rho_only(mp_type z, mp_type eta);
mp_type mu_only(mp_type z, mp_type eta);
mp_type rho_mu(mp_type z, mp_type eta);
mp_type rhop(mp_type z, mp_type eta);
mp_type psim(mp_type z, mp_type eta);
mp_type psip(mp_type z, mp_type eta);
mp_type rho_ex1(mp_type M, mp_type eta);
mp_type mu_example1(mp_type M, mp_type eta);
mp_type eta_eq(mp_type z, mp_type eta); // for example 1. need to find zero set of this for eta,z


// Note: values are turning complex here. Might be useful to use complex type. (Or its just fucked!)
mp_type C4_invariant(mp_type z, mp_type eta);

// Time Parameterization function
mp_type eta_to_time(mp_type z, mp_type eta);
mp_type et_function(mp_type z, mp_type eta);
mp_type ct_function(mp_type z, mp_type eta);
mp_type time_example1(mp_type M, mp_type eta);