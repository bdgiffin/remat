#include "gtest/gtest.h"
#include "ViscoElasticity.h"
#include "Parameters.h"
#include "types.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace {

using FloatViscoModel2D = ViscoElasticity<Real,Real>;

struct GradResult2D {
  Real loss;
  Real dparam_relaxation_time;
  Real dparam_shear_modulus_maxwell_element;
};

Parameters make_params_2d(Real tau, Real mu_e, int mat_overflow_limit = 1000000) {
  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 2.0;
  params["poissons_ratio"] = 0.25;
  params["relaxation_time"] = tau;
  params["shear_modulus_Maxwell_element"] = mu_e;
  params["mat_overflow_limit"] = mat_overflow_limit;
  return params;
}

void make_F_from_small_strain(const std::array<Real,3>& eps, Real (&F)[2][2]) {
  const Real exx = eps[0];
  const Real eyy = eps[1];
  const Real gxy = eps[2];
  F[0][0] = 1.0 + exx;
  F[1][1] = 1.0 + eyy;
  F[0][1] = 0.5*gxy;
  F[1][0] = 0.5*gxy;
}

GradResult2D run_objective_and_gradients_2d(Real tau, Real mu_e,
                                            const std::vector<std::array<Real,3> >& strain_history,
                                            Real dt, bool run_backward_adjoint) {
  Parameters params = make_params_2d(tau,mu_e,1000000);
  FloatViscoModel2D model(params);
  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  Real loss = 0.0;
  Real psi = 0.0;
  for (int i=0; i<int(strain_history.size()); i++) {
    const Real sxx = state[0];
    const Real syy = state[1];
    const Real sxy = state[5];
    loss += 0.5*(sxx*sxx + syy*syy + 2.0*sxy*sxy);

    Real F[2][2];
    make_F_from_small_strain(strain_history[i],F);
    model.update(F,psi,state.data(),dt,PassPhase::Forward);
  }

  if (run_backward_adjoint) {
    for (int i=int(strain_history.size())-1; i>=0; i--) {
      Real F[2][2];
      make_F_from_small_strain(strain_history[i],F);
      model.update(F,psi,state.data(),dt,PassPhase::BackwardAdjoint);
    }
  }

  return { loss, state[17], state[18] };
}

} // anonymous namespace

TEST(test_ViscoElasticityAdjoint, reversibility_primal_state_fields) {
  Parameters params = make_params_2d(0.3,0.8,1000000);
  FloatViscoModel2D model(params);

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());
  std::vector<Real> initial_state = state;

  std::vector<std::array<Real,3> > strain_history = {
    { 0.000, 0.000, 0.000 },
    { 0.010, 0.003, 0.006 },
    { 0.014, 0.005, 0.009 },
    { 0.009, 0.004, 0.004 },
    { 0.004, 0.001, 0.000 },
    { 0.000, 0.000, 0.000 }
  };

  const Real dt = 1.0e-2;
  Real psi = 0.0;
  for (int i=0; i<int(strain_history.size()); i++) {
    Real F[2][2];
    make_F_from_small_strain(strain_history[i],F);
    model.update(F,psi,state.data(),dt,PassPhase::Forward);
  }
  for (int i=int(strain_history.size())-1; i>=0; i--) {
    Real F[2][2];
    make_F_from_small_strain(strain_history[i],F);
    model.update(F,psi,state.data(),dt,PassPhase::Backward);
  }

  const Real tol = 2.0e-4;
  for (int i=0; i<=16; i++) {
    ASSERT_NEAR(state[i],initial_state[i],tol);
  }
}

TEST(test_ViscoElasticityAdjoint, backward_does_not_accumulate_gradients) {
  Parameters params = make_params_2d(0.3,0.8,1000000);
  FloatViscoModel2D model(params);

  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  std::vector<std::array<Real,3> > strain_history = {
    { 0.000, 0.000, 0.000 },
    { 0.010, 0.003, 0.006 },
    { 0.014, 0.005, 0.009 },
    { 0.009, 0.004, 0.004 }
  };

  const Real dt = 1.0e-2;
  Real psi = 0.0;
  for (int i=0; i<int(strain_history.size()); i++) {
    Real F[2][2];
    make_F_from_small_strain(strain_history[i],F);
    model.update(F,psi,state.data(),dt,PassPhase::Forward);
  }

  ASSERT_NEAR(state[17],0.0,1.0e-14);
  ASSERT_NEAR(state[18],0.0,1.0e-14);

  for (int i=int(strain_history.size())-1; i>=0; i--) {
    Real F[2][2];
    make_F_from_small_strain(strain_history[i],F);
    model.update(F,psi,state.data(),dt,PassPhase::Backward);
  }

  ASSERT_NEAR(state[17],0.0,1.0e-14);
  ASSERT_NEAR(state[18],0.0,1.0e-14);
}

TEST(test_ViscoElasticityAdjoint, gradients_match_finite_difference) {
  const Real tau = 0.35;
  const Real mu_e = 0.9;
  const Real dt = 2.0e-3;
  const Real h_tau = 5.0e-5;
  const Real h_mu = 5.0e-5;

  std::vector<std::array<Real,3> > strain_history = {
    { 0.000, 0.000, 0.000 },
    { 0.004, 0.001, 0.002 },
    { 0.008, 0.003, 0.004 },
    { 0.012, 0.005, 0.006 },
    { 0.010, 0.004, 0.004 },
    { 0.007, 0.002, 0.001 },
    { 0.003, 0.001, 0.000 },
    { 0.001, 0.000, 0.000 }
  };

  GradResult2D adj = run_objective_and_gradients_2d(tau,mu_e,strain_history,dt,true);
  GradResult2D plus_tau = run_objective_and_gradients_2d(tau + h_tau,mu_e,strain_history,dt,false);
  GradResult2D minus_tau = run_objective_and_gradients_2d(tau - h_tau,mu_e,strain_history,dt,false);
  Real fd_tau = (plus_tau.loss - minus_tau.loss)/(2.0*h_tau);

  Real tol_tau = 5.0e-2 * std::max(Real(1.0),std::fabs(fd_tau));
  ASSERT_NEAR(adj.dparam_relaxation_time,fd_tau,tol_tau);

  GradResult2D plus_mu = run_objective_and_gradients_2d(tau,mu_e + h_mu,strain_history,dt,false);
  GradResult2D minus_mu = run_objective_and_gradients_2d(tau,mu_e - h_mu,strain_history,dt,false);
  Real fd_mu = (plus_mu.loss - minus_mu.loss)/(2.0*h_mu);

  Real tol_mu = 5.0e-2 * std::max(Real(1.0),std::fabs(fd_mu));
  ASSERT_NEAR(adj.dparam_shear_modulus_maxwell_element,fd_mu,tol_mu);
}
