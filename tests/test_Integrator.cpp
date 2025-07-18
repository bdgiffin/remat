#include "gtest/gtest.h"
#include "Integrator.h"
#include "Rational.h"
#include "Fixed.h"
#include "Dual.h"
#include <vector>
#include <iostream>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_V = -6;
const int     EXPONENT_U = -4;
typedef Fixed<RADIX,EXPONENT_V> FixedV;
typedef Fixed<RADIX,EXPONENT_U> FixedU;

TEST(test_Integrator, constructors) {
  // Test creation of integrator
  Integrator<Rational> leapfrog;
} /* TEST(test_Integrator, constructors) */

TEST(test_Integrator, half_step_velocity_update) {

  const int N = 10;
  Real dt = 0.1;

  //// Test integrator on Real types
  //{
  //  Integrator<Real> leapfrog;
  //  std::vector<Real> v(N,0.0);
  //  std::vector<Real> u(N,0.0);
  //  std::vector<Real> a(N,0.0);
  //  std::vector<Real> alpha(N,1.0);
  //
  //  leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
  //  leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  //}
  //
  //// Test integrator on primal Fixed types
  //{
  //  Integrator<Rational> leapfrog;
  //  std::vector<FixedV> v(N,0.0);
  //  std::vector<FixedU> u(N,0.0);
  //  std::vector<Real> a(N,0.0);
  //  std::vector<Real> alpha(N,1.0);
  //
  //  // Update should be reversible provided alpha is unity
  //  leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
  //  leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  //}

  // Test integrator on Dual Fixed types
  {
    Integrator<Rational> leapfrog;
    std::vector<Dual<FixedV> > v(N,Dual<FixedV>(0.0,0.0));
    std::vector<Dual<FixedU> > u(N,Dual<FixedU>(0.0,0.0));
    std::vector<Dual<Real> >   a(N,Dual<Real>(0.0,0.0));
    std::vector<Real> alpha(N,1.0);

    // Update should be reversible provided no overflow occurs
    leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
    leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  }
  
} /* TEST(test_Integrator, constructors) */
