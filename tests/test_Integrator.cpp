#include "gtest/gtest.h"
#include "Integrator.h"
#include "Rational.h"
#include "Fixed.h"
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
  Integrator<FixedV,FixedU,Rational> leapfrog;
  
  const int N = 10;
  Matrix<N> A;
  Matrix<N> B(2.0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      ASSERT_EQ(A(i,j), 0.0);
      ASSERT_EQ(B(i,j), 2.0);
    }
  }
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "B = " << std::endl << B << std::endl;
} /* TEST(test_Integrator, constructors) */

TEST(test_Integrator, half_step_velocity_update) {

  const int N = 10;
  Real dt = 0.1;
  std::vector<Real> a(N,0.0);
  std::vector<Real> alpha(N,1.0);

  // Test integrator on Real types
  {
    Integrator<Real,Real,Real> leapfrog;
    std::vector<Real> v(N,0.0);
    std::vector<Real> u(N,0.0);

    leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
    leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  }

  // Test integrator on primal Fixed types
  {
    Integrator<FixedV,FixedU,Rational> leapfrog;
    std::vector<FixedV> v(N,0.0);
    std::vector<FixedU> u(N,0.0);

    // Update should be reversible provided alpha is unity
    leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
    leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  }

  // Test integrator on Dual Fixed types
  {
    Integrator<FixedV,FixedU,Rational> leapfrog;
    std::vector<Dual<FixedV> > v(N,0.0);
    std::vector<Dual<FixedU> > u(N,0.0);

    // Update should be reversible provided no overflow occurs
    leapfrog.first_half_step_velocity_update( +dt,v.data(),a.data(),alpha.data(),N); // forward
    leapfrog.second_half_step_velocity_update(-dt,v.data(),a.data(),alpha.data(),N); // reverse
  }
  
  const int N = 10;
  Matrix<N> A;
  Matrix<N> B(2.0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      ASSERT_EQ(A(i,j), 0.0);
      ASSERT_EQ(B(i,j), 2.0);
    }
  }
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "B = " << std::endl << B << std::endl;
} /* TEST(test_Integrator, constructors) */
