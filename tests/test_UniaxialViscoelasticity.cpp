#include "gtest/gtest.h"
#include "UniaxialViscoelasticity.h"
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

using FloatViscoModel = UniaxialViscoelasticity<Real,Real,Real>;
using FixedViscoModel = UniaxialViscoelasticity<FixedE,Rational,Real>;

struct GradResult {
  Real loss;
  Real df_dtau;
  Real df_dE;
};

template<class ModelT>
Parameters make_params(Real tau, Real E, int mat_overflow_limit) {
  Parameters params;
  params["truss_density"] = 1.0;
  params["area"] = 1.0;
  params["truss_youngs_modulus"] = E;
  params["relaxation_time"] = tau;
  params["mat_overflow_limit"] = mat_overflow_limit;
  return params;
}

GradResult run_objective_and_gradients(Real tau, Real E,
				       const std::vector<Real>& strain_history,
				       Real dt, bool run_backward) {
  Parameters params = make_params<FloatViscoModel>(tau,E,1000000);
  FloatViscoModel model(params);
  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  Real loss = 0.0;
  Real psi = 0.0;
  for (int i = 0; i < int(strain_history.size()); ++i) {
    // The built-in adjoint recurrence uses sigma_n before the n->n+1 constitutive update.
    Real sigma_n = state[0];
    loss += 0.5 * sigma_n * sigma_n / E;

    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt);
  }

  if (run_backward) {
    for (int i = int(strain_history.size())-1; i >= 0; --i) {
      Real lambda = 1.0 + strain_history[i];
      model.update(lambda,psi,state.data(),-dt);
    }
  }

  return { loss, state[5], state[6] };
}

} // end anonymous namespace

TEST(test_UniaxialViscoelasticity, constructors) {
  // Empty constructor
  {
    FixedViscoModel model;
  }

  // Constructor from parameters
  {
    Parameters params = make_params<FixedViscoModel>(0.3,2.5,30);
    FixedViscoModel model(params);
  }
}

TEST(test_UniaxialViscoelasticity, reversibility_primal_state_fields) {
  Parameters params = make_params<FloatViscoModel>(0.25,2.0,1000000);
  FloatViscoModel model(params);

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());
  std::vector<Real> initial_state = state;

  std::vector<Real> strain_history = {
    0.0000, 0.0100, 0.0200, 0.0150, 0.0050, 0.0000, -0.0040, 0.0000
  };

  const Real dt = 1.0e-2;
  Real psi = 0.0;
  for (int i = 0; i < int(strain_history.size()); ++i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt);
  }
  for (int i = int(strain_history.size())-1; i >= 0; --i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),-dt);
  }

  // Primal and reversible bookkeeping fields should return to the initial state.
  const Real tol = 2.0e-4;
  ASSERT_NEAR(state[0],initial_state[0],tol); // axial_stress
  ASSERT_NEAR(state[1],initial_state[1],tol); // axial_strain
  ASSERT_NEAR(state[2],initial_state[2],tol); // viscous_strain
  ASSERT_NEAR(state[3],initial_state[3],tol); // dual_viscous_strain
  ASSERT_NEAR(state[4],initial_state[4],tol); // overflow_counter
}

TEST(test_UniaxialViscoelasticity, gradient_df_dtau_matches_finite_difference) {
  const Real tau = 0.35;
  const Real E   = 2.0;
  const Real dt  = 2.5e-3;
  const Real h_tau = 1.0e-4;

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0100, 0.0140, 0.0120, 0.0080,
    0.0040, 0.0010, -0.0010, 0.0000, 0.0020, 0.0040
  };

  GradResult adj = run_objective_and_gradients(tau,E,strain_history,dt,true);
  GradResult plus  = run_objective_and_gradients(tau + h_tau,E,strain_history,dt,false);
  GradResult minus = run_objective_and_gradients(tau - h_tau,E,strain_history,dt,false);
  Real fd_dtau = (plus.loss - minus.loss)/(2.0*h_tau);

  Real tol = 1.0e-2 * std::max(Real(1.0),std::fabs(fd_dtau));
  ASSERT_NEAR(adj.df_dtau,fd_dtau,tol);
}

TEST(test_UniaxialViscoelasticity, gradient_df_dE_matches_finite_difference) {
  const Real tau = 0.35;
  const Real E   = 2.0;
  const Real dt  = 2.5e-3;
  const Real h_E = 1.0e-4;

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0100, 0.0140, 0.0120, 0.0080,
    0.0040, 0.0010, -0.0010, 0.0000, 0.0020, 0.0040
  };

  GradResult adj = run_objective_and_gradients(tau,E,strain_history,dt,true);
  GradResult plus  = run_objective_and_gradients(tau,E + h_E,strain_history,dt,false);
  GradResult minus = run_objective_and_gradients(tau,E - h_E,strain_history,dt,false);
  Real fd_dE = (plus.loss - minus.loss)/(2.0*h_E);

  Real tol = 1.0e-2 * std::max(Real(1.0),std::fabs(fd_dE));
  ASSERT_NEAR(adj.df_dE,fd_dE,tol);
}

TEST(test_UniaxialViscoelasticity, explicit_mode_forward_matches_legacy_forward) {
  Parameters params = make_params<FloatViscoModel>(0.35,2.0,1000000);
  FloatViscoModel model(params);

  std::vector<Real> state_legacy(model.num_state_vars(),0.0);
  std::vector<Real> state_mode(model.num_state_vars(),0.0);
  model.initialize(state_legacy.data());
  model.initialize(state_mode.data());

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0100, 0.0140, 0.0120, 0.0080,
    0.0040, 0.0010, -0.0010, 0.0000, 0.0020, 0.0040
  };

  const Real dt = 2.5e-3;
  Real psi_legacy = 0.0;
  Real psi_mode = 0.0;

  for (int i = 0; i < int(strain_history.size()); ++i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi_legacy,state_legacy.data(),dt);
    model.update(lambda,psi_mode,state_mode.data(),dt,FloatViscoModel::MaterialUpdateMode::Forward);
  }

  const Real tol = 1.0e-12;
  for (int i = 0; i < model.num_state_vars(); ++i) {
    ASSERT_NEAR(state_legacy[i],state_mode[i],tol);
  }
  ASSERT_NEAR(psi_legacy,psi_mode,tol);
}

TEST(test_UniaxialViscoelasticity, explicit_remat_and_adjoint_match_legacy_backward) {
  Parameters params = make_params<FloatViscoModel>(0.35,2.0,1000000);
  FloatViscoModel model(params);

  std::vector<Real> state_legacy(model.num_state_vars(),0.0);
  std::vector<Real> state_split(model.num_state_vars(),0.0);
  model.initialize(state_legacy.data());
  model.initialize(state_split.data());

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0100, 0.0140, 0.0120, 0.0080,
    0.0040, 0.0010, -0.0010, 0.0000, 0.0020, 0.0040
  };

  const Real dt = 2.5e-3;
  Real psi = 0.0;
  for (int i = 0; i < int(strain_history.size()); ++i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state_legacy.data(),dt);
    model.update(lambda,psi,state_split.data(),dt);
  }

  for (int i = int(strain_history.size())-1; i >= 0; --i) {
    Real lambda = 1.0 + strain_history[i];

    // Legacy combined backward/remat + adjoint.
    model.update(lambda,psi,state_legacy.data(),-dt);

    // Explicit split: local adjoint pass first (without changing primal state),
    // then rematerialize primal state.
    model.update(lambda,psi,state_split.data(),dt,FloatViscoModel::MaterialUpdateMode::AdjointBackward);
    model.update(lambda,psi,state_split.data(),dt,FloatViscoModel::MaterialUpdateMode::RematBackward);
  }

  const Real tol = 1.0e-10;
  for (int i = 0; i < model.num_state_vars(); ++i) {
    ASSERT_NEAR(state_legacy[i],state_split[i],tol);
  }
}

TEST(test_UniaxialViscoelasticity, external_seed_matches_default_objective_policy) {
  Parameters params = make_params<FloatViscoModel>(0.35,2.0,1000000);
  FloatViscoModel model(params);
  const Real E = params["truss_youngs_modulus"];

  std::vector<Real> state_default(model.num_state_vars(),0.0);
  std::vector<Real> state_external(model.num_state_vars(),0.0);
  model.initialize(state_default.data());
  model.initialize(state_external.data());

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0100, 0.0140, 0.0120
  };

  const Real dt = 2.5e-3;
  Real psi = 0.0;
  for (int i = 0; i < int(strain_history.size()); ++i) {
    const Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state_default.data(),dt);
    model.update(lambda,psi,state_external.data(),dt);
  }

  // Legacy backward path with built-in objective seed.
  const int last = int(strain_history.size()) - 1;
  const Real lambda_last = 1.0 + strain_history[last];
  model.update(lambda_last,psi,state_default.data(),-dt);

  // Explicit adjoint path with external objective seed.
  std::vector<Real> state_n_tmp = state_external;
  model.update(lambda_last,psi,state_n_tmp.data(),dt,FloatViscoModel::MaterialUpdateMode::RematBackward);
  const Real sigma_n = state_n_tmp[0];

  FloatViscoModel::LocalAdjointSeed seed;
  seed.bar_sigma_n = sigma_n/E;
  seed.bar_sigma_np1 = 0.0;
  seed.direct_dE = -0.5*(sigma_n/E)*(sigma_n/E);
  model.set_external_objective_seed(seed);

  model.update(lambda_last,psi,state_external.data(),dt,FloatViscoModel::MaterialUpdateMode::AdjointBackward);
  model.update(lambda_last,psi,state_external.data(),dt,FloatViscoModel::MaterialUpdateMode::RematBackward);

  const Real tol = 1.0e-12;
  ASSERT_NEAR(state_default[5],state_external[5],tol); // df_dtau
  ASSERT_NEAR(state_default[6],state_external[6],tol); // df_dE
  ASSERT_NEAR(state_default[7],state_external[7],tol); // lambda_adjoint
}
