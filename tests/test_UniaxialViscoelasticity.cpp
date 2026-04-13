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
    model.update(lambda,psi,state.data(),dt,PassPhase::Forward);
  }

  if (run_backward) {
    for (int i = int(strain_history.size())-1; i >= 0; --i) {
      Real lambda = 1.0 + strain_history[i];
      model.update(lambda,psi,state.data(),dt,PassPhase::BackwardAdjoint);
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
    model.update(lambda,psi,state.data(),dt,PassPhase::Forward);
  }
  for (int i = int(strain_history.size())-1; i >= 0; --i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt,PassPhase::Backward);
  }

  // Primal and reversible bookkeeping fields should return to the initial state.
  const Real tol = 2.0e-4;
  ASSERT_NEAR(state[0],initial_state[0],tol); // axial_stress
  ASSERT_NEAR(state[1],initial_state[1],tol); // axial_strain
  ASSERT_NEAR(state[2],initial_state[2],tol); // viscous_strain
  ASSERT_NEAR(state[3],initial_state[3],tol); // dual_viscous_strain
  ASSERT_NEAR(state[4],initial_state[4],tol); // overflow_counter
}

TEST(test_UniaxialViscoelasticity, backward_does_not_accumulate_gradients) {
  Parameters params = make_params<FloatViscoModel>(0.25,2.0,1000000);
  FloatViscoModel model(params);

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  std::vector<Real> strain_history = {
    0.0000, 0.0100, 0.0200, 0.0150, 0.0050, 0.0000
  };

  const Real dt = 1.0e-2;
  Real psi = 0.0;
  for (int i = 0; i < int(strain_history.size()); ++i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt,PassPhase::Forward);
  }

  ASSERT_NEAR(state[5],0.0,1.0e-14);
  ASSERT_NEAR(state[6],0.0,1.0e-14);

  for (int i = int(strain_history.size())-1; i >= 0; --i) {
    Real lambda = 1.0 + strain_history[i];
    model.update(lambda,psi,state.data(),dt,PassPhase::Backward);
  }

  ASSERT_NEAR(state[5],0.0,1.0e-14);
  ASSERT_NEAR(state[6],0.0,1.0e-14);
}

TEST(test_UniaxialViscoelasticity, backward_and_backwardadjoint_match_primal_trajectory) {
  Parameters params = make_params<FloatViscoModel>(0.35,2.3,1000000);
  FloatViscoModel model_backward(params);
  FloatViscoModel model_adjoint(params);

  std::vector<Real> state_backward(model_backward.num_state_vars(),0.0);
  std::vector<Real> state_adjoint(model_adjoint.num_state_vars(),0.0);
  model_backward.initialize(state_backward.data());
  model_adjoint.initialize(state_adjoint.data());

  std::vector<Real> strain_history = {
    0.0000, 0.0050, 0.0120, 0.0200, 0.0170, 0.0100, 0.0020, -0.0010
  };

  const Real dt = 1.0e-2;
  Real psi_backward = 0.0;
  Real psi_adjoint = 0.0;

  for (int i = 0; i < int(strain_history.size()); ++i) {
    Real lambda = 1.0 + strain_history[i];
    model_backward.update(lambda,psi_backward,state_backward.data(),dt,PassPhase::Forward);
    model_adjoint.update(lambda,psi_adjoint,state_adjoint.data(),dt,PassPhase::Forward);
  }

  const Real tol = 2.0e-4;
  for (int i = int(strain_history.size())-1; i >= 0; --i) {
    Real lambda = 1.0 + strain_history[i];
    model_backward.update(lambda,psi_backward,state_backward.data(),dt,PassPhase::Backward);
    model_adjoint.update(lambda,psi_adjoint,state_adjoint.data(),dt,PassPhase::BackwardAdjoint);

    ASSERT_NEAR(state_backward[0],state_adjoint[0],tol); // axial_stress
    ASSERT_NEAR(state_backward[1],state_adjoint[1],tol); // axial_strain
    ASSERT_NEAR(state_backward[2],state_adjoint[2],tol); // viscous_strain
    ASSERT_NEAR(state_backward[3],state_adjoint[3],tol); // dual_viscous_strain
    ASSERT_NEAR(state_backward[4],state_adjoint[4],tol); // overflow_counter
  }
}

TEST(test_UniaxialViscoelasticity, local_forward_backward_closed_form) {
  const Real tau = 0.4;
  const Real E = 2.5;
  const Real dt = 1.0e-2;
  Parameters params = make_params<FloatViscoModel>(tau,E,1000000);
  FloatViscoModel model(params);

  FloatViscoModel::LocalForwardInput in = { 0.08, 0.03, dt };
  FloatViscoModel::LocalForwardOutput out = model.forward_local(in);

  const Real A = std::exp(-dt/tau);
  const Real q_expected = A*in.q_n + (1.0-A)*in.eps_np1;
  const Real sigma_expected = E*(in.eps_np1 - q_expected);
  ASSERT_NEAR(out.A,A,1.0e-14);
  ASSERT_NEAR(out.q_np1,q_expected,1.0e-14);
  ASSERT_NEAR(out.sigma_np1,sigma_expected,1.0e-14);

  FloatViscoModel::LocalBackwardInput bin;
  bin.in = in;
  bin.out = out;
  bin.bar_sigma_np1 = 0.7;
  bin.bar_q_np1 = -0.25;
  FloatViscoModel::LocalBackwardOutput bout = model.backward_local(bin);

  const Real bar_eps_expected = E*A*bin.bar_sigma_np1 + (1.0-A)*bin.bar_q_np1;
  const Real bar_qn_expected = -E*A*bin.bar_sigma_np1 + A*bin.bar_q_np1;
  const Real dE_expected = A*(in.eps_np1 - in.q_n)*bin.bar_sigma_np1;
  const Real dTau_expected =
    (in.eps_np1 - in.q_n)*(E*bin.bar_sigma_np1 - bin.bar_q_np1)*A*(dt/(tau*tau));

  ASSERT_NEAR(bout.bar_eps_np1,bar_eps_expected,1.0e-14);
  ASSERT_NEAR(bout.bar_q_n,bar_qn_expected,1.0e-14);
  ASSERT_NEAR(bout.grad.dE(),dE_expected,1.0e-14);
  ASSERT_NEAR(bout.grad.dTau(),dTau_expected,1.0e-14);
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

  const Real h_E = 1.0e-4;
  GradResult plus_E  = run_objective_and_gradients(tau,E + h_E,strain_history,dt,false);
  GradResult minus_E = run_objective_and_gradients(tau,E - h_E,strain_history,dt,false);
  Real fd_dE = (plus_E.loss - minus_E.loss)/(2.0*h_E);

  Real tol_E = 1.0e-2 * std::max(Real(1.0),std::fabs(fd_dE));
  ASSERT_NEAR(adj.df_dE,fd_dE,tol_E);
}
