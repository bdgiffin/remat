#include "gtest/gtest.h"
#include "UniaxialViscoplasticity.h"
#include "Rational.h"
#include "Fixed.h"
#include "Parameters.h"
#include "types.h"
#include <vector>
#include <cmath>
#include <algorithm>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_E = -6;
typedef Fixed<RADIX,EXPONENT_E> FixedE;

namespace {

using FloatViscoplasticModel = UniaxialViscoplasticity<Real,Real,Real>;
using FixedViscoplasticModel = UniaxialViscoplasticity<FixedE,Rational,Real>;

struct GradResult {
  Real loss = 0.0;
  Real dL_drelaxation_time = 0.0;
};

Parameters make_params(Real tau, Real E, Real yield_stress, int mat_overflow_limit) {
  Parameters params;
  params["truss_density"] = 1.0;
  params["area"] = 1.0;
  params["truss_youngs_modulus"] = E;
  params["relaxation_time"] = tau;
  params["yield_stress"] = yield_stress;
  params["eps_fail"] = 10.0;
  params["mat_overflow_limit"] = mat_overflow_limit;
  return params;
}

GradResult run_objective_and_gradients(Real tau, Real E, Real yield_stress,
                                       const std::vector<Real>& strain_history,
                                       Real dt, bool run_backward) {
  Parameters params = make_params(tau,E,yield_stress,1000000);
  FloatViscoplasticModel model(params);
  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  Real loss = 0.0;
  Real psi = 0.0;
  for (int i=0; i<int(strain_history.size()); ++i) {
    const Real sigma_n = model.axial_stress(state.data());
    loss += 0.5*sigma_n*sigma_n/E;
    const Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt);
  }

  if (run_backward) {
    for (int i=int(strain_history.size())-1; i>=0; --i) {
      const Real lambda = 1.0 + strain_history[i];
      model.update(lambda,psi,state.data(),-dt);
    }
  }

  return {loss,state[7]};
}

} // namespace

TEST(test_UniaxialViscoplasticity, constructors) {
  {
    FixedViscoplasticModel model;
  }
  {
    Parameters params = make_params(0.2,2000.0,30.0,50);
    FixedViscoplasticModel model(params);
  }
}

TEST(test_UniaxialViscoplasticity, reversibility_primal_state_fields) {
  Parameters params = make_params(0.3,2000.0,30.0,1000000);
  FloatViscoplasticModel model(params);

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());
  std::vector<Real> initial_state = state;

  const std::vector<Real> strain_history = {
    0.0000, 0.0120, 0.0200, 0.0280, 0.0350, 0.0250, 0.0100, 0.0000,
   -0.0100,-0.0200,-0.0280,-0.0200,-0.0100, 0.0000
  };

  const Real dt = 7.0e-3;
  Real psi = 0.0;
  for (int i=0; i<int(strain_history.size()); ++i) {
    const Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt);
  }
  for (int i=int(strain_history.size())-1; i>=0; --i) {
    const Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),-dt);
  }

  const Real tol = 3.0e-4;
  for (int i=0; i<=6; ++i) {
    ASSERT_NEAR(state[i],initial_state[i],tol);
  }
}

TEST(test_UniaxialViscoplasticity, gradient_dL_drelaxation_time_matches_finite_difference) {
  const Real tau = 0.35;
  const Real E = 2.0;
  const Real yield_stress = 0.02;
  const Real dt = 2.0e-3;
  const Real h_tau = 2.0e-4;

  const std::vector<Real> strain_history = {
    0.0000, 0.0120, 0.0150, 0.0180, 0.0210, 0.0240,
    0.0270, 0.0300, 0.0320, 0.0290, 0.0260, 0.0230
  };

  const GradResult adj = run_objective_and_gradients(tau,E,yield_stress,strain_history,dt,true);
  const GradResult plus = run_objective_and_gradients(tau+h_tau,E,yield_stress,strain_history,dt,false);
  const GradResult minus = run_objective_and_gradients(tau-h_tau,E,yield_stress,strain_history,dt,false);
  const Real fd_dtau = (plus.loss - minus.loss)/(2.0*h_tau);

  const Real tol = 6.0e-2*std::max(Real(1.0),std::fabs(fd_dtau));
  ASSERT_NEAR(adj.dL_drelaxation_time,fd_dtau,tol);
}
