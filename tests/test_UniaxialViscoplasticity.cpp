#include "gtest/gtest.h"
#include "UniaxialViscoplasticity.h"
#include "Rational.h"
#include "Fixed.h"
#include "Dual.h"
#include "Parameters.h"
#include "types.h"
#include <vector>
#include <iostream>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_E = -6;
typedef Fixed<RADIX,EXPONENT_E> FixedE;

TEST(test_UniaxialViscoplasticity, constructors) {
  // Test creation of UniaxialViscoplasticity model
  {
    UniaxialViscoplasticity<FixedE,Rational> model;
  }

  // Test creation of UniaxialViscoplasticity model from parameters
  {
    Parameters params;
    params["truss_density"] = 1.0;
    params["area"] = 1.0;
    params["truss_youngs_modulus"] = 2000.0;
    params["viscosity"] = 1.0e+1;
    params["yield_stress"] = 30.0;
    params["eps_fail"] = 0.3;
    params["mat_overflow_limit"] = 50;
    UniaxialViscoplasticity<FixedE,Rational> model(params);
  }
} /* TEST(test_UniaxialViscoplasticity, constructors) */

TEST(test_UniaxialViscoplasticity, update) {
  // Initialize model
  Parameters params;
  params["truss_density"] = 1.0;
  params["area"] = 1.0;
  params["truss_youngs_modulus"] = 2000.0;
  params["viscosity"] = 1.0e+1;
  params["yield_stress"] = 30.0;
  params["eps_fail"] = 0.3;
  params["mat_overflow_limit"] = 50;
  FixedE yield_strain(params["yield_stress"]/params["truss_youngs_modulus"]);
  UniaxialViscoplasticity<FixedE,Rational> model(params);
  
  Real deps = 1.0e-6; // strain increment
  Real dt   = 7.0e-3; // time step
  Real psi;           // internal energy density
  int Nsteps = 10;

  // Loading/unloading in tension
  {

    // Initialize model state
    std::vector<Real> state(model.num_state_vars());
    model.initialize(state.data());

    // Update model state
    Real lambda = 1.0 + Real(yield_strain) - deps; // stretch ratio
    Real time = 0.0;
    std::vector<Real> initial_state;
    for (int i = 0; i <= Nsteps; i++) {
      model.update(lambda, psi, state.data(), dt);
      if (i == 0) initial_state = state;
      FixedE plastic_strain, dual_plastic_strain;
      load_from_Real(state[2],plastic_strain);
      load_from_Real(state[3],dual_plastic_strain);
      std::cout << "Step " << i << std::endl;
      std::cout << "                 time = " << time << std::endl;
      std::cout << "         total strain = " << lambda - 1.0 << std::endl;
      std::cout << "       plastic strain = " << plastic_strain << std::endl;
      std::cout << "  dual plastic strain = " << dual_plastic_strain << std::endl;
      lambda += deps*(1.0 - 2*(i*i)/Real(Nsteps*Nsteps));
      time   += dt;
    }

    // Reverse update model state
    lambda -= deps*(1.0 - 2*(Nsteps*Nsteps)/Real(Nsteps*Nsteps));
    time   -= dt;
    std::vector<Real> final_state;
    for (int i = Nsteps-1; i >= 0; i--) {
      lambda -= deps*(1.0 - 2*(i*i)/Real(Nsteps*Nsteps));
      time   -= dt;
      model.update(lambda, psi, state.data(), -dt);
      if (i == 0) final_state = state;
      FixedE plastic_strain, dual_plastic_strain;
      load_from_Real(state[2],plastic_strain);
      load_from_Real(state[3],dual_plastic_strain);
      std::cout << "Step " << i << std::endl;
      std::cout << "                 time = " << time << std::endl;
      std::cout << "         total strain = " << lambda - 1.0 << std::endl;
      std::cout << "       plastic strain = " << plastic_strain << std::endl;
      std::cout << "  dual plastic strain = " << dual_plastic_strain << std::endl;
    }

    for (int i = 0; i < model.num_state_vars(); i++) {
      ASSERT_EQ(initial_state[i],final_state[i]);
    }
  }

  // Loading/unloading in compression
  {

    // Initialize model state
    std::vector<Real> state(model.num_state_vars());
    model.initialize(state.data());

    // Update model state
    Real lambda = 1.0 - Real(yield_strain) + deps; // stretch ratio
    Real time = 0.0;
    std::vector<Real> initial_state;
    for (int i = 0; i <= Nsteps; i++) {
      model.update(lambda, psi, state.data(), dt);
      if (i == 0) initial_state = state;
      FixedE plastic_strain, dual_plastic_strain;
      load_from_Real(state[2],plastic_strain);
      load_from_Real(state[3],dual_plastic_strain);
      std::cout << "Step " << i << std::endl;
      std::cout << "                 time = " << time << std::endl;
      std::cout << "         total strain = " << lambda - 1.0 << std::endl;
      std::cout << "       plastic strain = " << plastic_strain << std::endl;
      std::cout << "  dual plastic strain = " << dual_plastic_strain << std::endl;
      lambda -= deps*(1.0 - 2*(i*i)/Real(Nsteps*Nsteps));
      time   += dt;
    }

    // Reverse update model state
    lambda += deps*(1.0 - 2*(Nsteps*Nsteps)/Real(Nsteps*Nsteps));
    time   -= dt;
    std::vector<Real> final_state;
    for (int i = Nsteps-1; i >= 0; i--) {
      lambda += deps*(1.0 - 2*(i*i)/Real(Nsteps*Nsteps));
      time   -= dt;
      model.update(lambda, psi, state.data(), -dt);
      if (i == 0) final_state = state;
      FixedE plastic_strain, dual_plastic_strain;
      load_from_Real(state[2],plastic_strain);
      load_from_Real(state[3],dual_plastic_strain);
      std::cout << "Step " << i << std::endl;
      std::cout << "                 time = " << time << std::endl;
      std::cout << "         total strain = " << lambda - 1.0 << std::endl;
      std::cout << "       plastic strain = " << plastic_strain << std::endl;
      std::cout << "  dual plastic strain = " << dual_plastic_strain << std::endl;
    }

    for (int i = 0; i < model.num_state_vars(); i++) {
      ASSERT_EQ(initial_state[i],final_state[i]);
    }
  }
  
} /* TEST(test_UniaxialViscoplasticity, update) */

TEST(test_UniaxialViscoplasticity, adjoint_scaffold_reports_unsupported) {
  Parameters params;
  params["truss_density"] = 1.0;
  params["area"] = 1.0;
  params["truss_youngs_modulus"] = 2000.0;
  params["viscosity"] = 1.0e+1;
  params["yield_stress"] = 30.0;
  params["eps_fail"] = 0.3;
  params["mat_overflow_limit"] = 50;
  UniaxialViscoplasticity<FixedE,Rational> model(params);

  EXPECT_EQ(model.adjoint_support_level(),int(AdjointSupportLevel::ScaffoldOnly));
  EXPECT_STREQ(model.adjoint_support_status(),"scaffold_only");

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  ScalarStepInput input;
  input.strain_n = 0.0;
  input.strain_np1 = 1.0e-3;
  input.dt = 1.0e-3;

  ScalarLocalSeed seed;
  seed.bar_sigma_n = 1.0;
  seed.bar_sigma_np1 = 0.0;
  seed.direct_dE = 0.0;

  UniaxialViscoplasticity<FixedE,Rational>::GradientAccumulator grad;
  auto response = model.adjoint_update(input,state.data(),seed,grad);

  EXPECT_FALSE(response.supported);
  EXPECT_EQ(grad.df_dtau,0.0);
  EXPECT_EQ(grad.df_dE,0.0);
  EXPECT_FALSE(grad.supported);
}
