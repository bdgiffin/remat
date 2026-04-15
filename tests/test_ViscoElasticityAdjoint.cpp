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
  Real dparam_stiffness_scaling_factor;
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

  return {
    loss,
    model.adjoint_get_param_gradient(state.data(),0),
    model.adjoint_get_param_gradient(state.data(),1),
    model.adjoint_get_param_gradient(state.data(),2)
  };
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

  // Finite-difference check for distributed stiffness scaling factor at one material point.
  // Use a constant perturbation of the internal state stiffness-scaling factor (equivalent to dL/ds).
  const Real h_scale = 5.0e-5;
  auto run_with_scale_offset = [&](Real scale_offset) {
    Parameters params = make_params_2d(tau,mu_e,1000000);
    FloatViscoModel2D model(params);
    std::vector<Real> state(model.num_state_vars(),0.0);
    model.initialize(state.data());
    state[15] = 1.0 + scale_offset; // STIFFNESS_SCALING field

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
    return loss;
  };

  const Real loss_plus_scale = run_with_scale_offset(+h_scale);
  const Real loss_minus_scale = run_with_scale_offset(-h_scale);
  const Real fd_scale = (loss_plus_scale - loss_minus_scale)/(2.0*h_scale);
  const Real tol_scale = 7.0e-2 * std::max(Real(1.0),std::fabs(fd_scale));
  ASSERT_NEAR(adj.dparam_stiffness_scaling_factor,fd_scale,tol_scale);
}

TEST(test_ViscoElasticityAdjoint, single_step_stiffness_scaling_gradient_matches_finite_difference) {
  const Real tau = 0.25;
  const Real mu_e = 0.7;
  const Real dt = 1.0e-3;
  const Real h_scale = 5.0e-5;

  std::vector<std::array<Real,3> > strain_history = {
    { 0.006, 0.002, 0.003 }
  };

  Parameters params = make_params_2d(tau,mu_e,1000000);
  FloatViscoModel2D model(params);
  std::vector<Real> state(model.num_state_vars(),0.0);
  model.initialize(state.data());

  Real psi = 0.0;
  Real F[2][2];
  make_F_from_small_strain(strain_history[0],F);
  model.update(F,psi,state.data(),dt,PassPhase::Forward);
  model.update(F,psi,state.data(),dt,PassPhase::BackwardAdjoint);
  const Real adj_scale = model.adjoint_get_param_gradient(state.data(),2);

  auto one_step_loss = [&](Real scale_offset) {
    Parameters p = make_params_2d(tau,mu_e,1000000);
    FloatViscoModel2D m(p);
    std::vector<Real> st(m.num_state_vars(),0.0);
    m.initialize(st.data());
    st[15] = 1.0 + scale_offset; // STIFFNESS_SCALING
    Real local_psi = 0.0;
    Real local_F[2][2];
    make_F_from_small_strain(strain_history[0],local_F);
    m.update(local_F,local_psi,st.data(),dt,PassPhase::Forward);
    const Real sxx = st[0];
    const Real syy = st[1];
    const Real sxy = st[5];
    return 0.5*(sxx*sxx + syy*syy + 2.0*sxy*sxy);
  };

  const Real fd_scale = (one_step_loss(+h_scale) - one_step_loss(-h_scale))/(2.0*h_scale);
  const Real tol = 8.0e-2 * std::max(Real(1.0),std::fabs(fd_scale));
  ASSERT_NEAR(adj_scale,fd_scale,tol);
}

TEST(test_ViscoElasticityAdjoint, second_kick_history_seed_enables_tau_gradient_match) {
  const Real tau = 0.28;
  const Real mu_e = 0.9;
  const Real dt = 1.0e-3;
  const Real h_tau = 5.0e-5;
  const Real seed_xx = 0.7;
  const Real seed_yy = -0.5;
  const Real seed_xy = 0.4;
  const std::array<Real,3> eps = { 0.020, 0.006, 0.010 };

  auto objective_one_step = [&](Real tau_value) {
    Parameters params = make_params_2d(tau_value,mu_e,1000000);
    params["adjoint_material_objective_weight"] = 0.0;
    FloatViscoModel2D model(params);
    std::vector<Real> state(model.num_state_vars(),0.0);
    model.initialize(state.data());
    Real F[2][2];
    make_F_from_small_strain(eps,F);
    Real psi = 0.0;
    model.update(F,psi,state.data(),dt,PassPhase::Forward);
    return seed_xx*state[0] + seed_yy*state[1] + seed_xy*state[5];
  };

  auto run_adjoint_tau = [&](bool with_history_seed) {
    Parameters params = make_params_2d(tau,mu_e,1000000);
    params["adjoint_material_objective_weight"] = 0.0;
    FloatViscoModel2D model(params);
    std::vector<Real> state(model.num_state_vars(),0.0);
    model.initialize(state.data());
    Real F[2][2];
    make_F_from_small_strain(eps,F);
    Real psi = 0.0;
    model.update(F,psi,state.data(),dt,PassPhase::Forward);
    if (with_history_seed) {
      model.adjoint_add_history_seed_from_stress(state.data(),seed_xx,seed_yy,seed_xy);
    }
    model.update(F,psi,state.data(),dt,PassPhase::BackwardAdjoint);
    return model.adjoint_get_param_gradient(state.data(),0);
  };

  const Real fd_tau = (objective_one_step(tau + h_tau) - objective_one_step(tau - h_tau))/(2.0*h_tau);
  const Real grad_without_history = run_adjoint_tau(false);
  const Real grad_with_history = run_adjoint_tau(true);
  const Real rel_without = std::fabs(grad_without_history - fd_tau) /
    std::max({std::fabs(grad_without_history),std::fabs(fd_tau),Real(1.0e-14)});
  const Real rel_with = std::fabs(grad_with_history - fd_tau) /
    std::max({std::fabs(grad_with_history),std::fabs(fd_tau),Real(1.0e-14)});

  ASSERT_GT(rel_without,0.5);
  ASSERT_LT(rel_with,5.0e-2);
}
