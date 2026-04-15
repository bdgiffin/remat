#include "gtest/gtest.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "ViscoElasticity.h"
#include "Truss.h"
#include "UniaxialMaterial.h"
#include "UniaxialViscoelasticity.h"
#include "Parameters.h"
#include "types.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace {

using TrussSystem = System<
  Element<Material>,
  Truss<UniaxialViscoelasticity<Real,Real,Real> >,
  Real, Real, Real>;

using ElementSystem = System<
  Element<ViscoElasticity<Real,Real> >,
  Truss<UniaxialMaterial>,
  Real, Real, Real>;

Real global_field_by_name(SystemBase& sys, const std::string& name) {
  const int n = sys.get_num_fields("global");
  std::vector<double> values(n,0.0);
  sys.get_fields("global",values.data());
  for (int i=0; i<n; i++) {
    if (name == sys.get_field_name("global",i)) return values[i];
  }
  return 0.0;
}

struct TrussKickRun {
  Real grad_tau;
  Real grad_E;
  Real phi;
};

TrussKickRun run_truss_case(Real tau, Real E, int steps, Real dt, Real seed, bool run_adjoint) {
  TrussSystem sys;

  const int Nnodes = 2;
  const int Ndofs_per_node = 2;
  const int Nelems = 0;
  const int Nnodes_per_elem = 4;

  double coordinates[4] = { 0.0,0.0, 1.0,0.0 };
  double velocities[4]  = { 0.0,0.0, 0.12,0.0 };
  bool fixity[4]        = { true,true, false,true };
  int connectivity_dummy[4] = { 0,0,0,0 };
  int truss_connectivity[2] = { 0,1 };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;
  params["truss_density"] = 1.0;
  params["area"] = 1.0;
  params["truss_youngs_modulus"] = E;
  params["relaxation_time"] = tau;
  params["mass_damping_factor"] = 0.0;
  params["contact_stiffness"] = 0.0;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity_dummy,Nelems,Nnodes_per_elem,params);
  sys.initialize_truss_elements(truss_connectivity,1,params);
  sys.initialize_state();

  for (int i=0; i<steps; i++) {
    sys.update_state(dt,PassPhase::Forward);
  }

  const int seed_dof = 2; // node 1, x
  Real phi = seed*Real(sys.v[seed_dof].first);

  Real grad_tau = 0.0;
  Real grad_E = 0.0;
  if (run_adjoint) {
    std::fill(sys.u_adjoint.begin(),sys.u_adjoint.end(),0.0);
    std::fill(sys.v_adjoint.begin(),sys.v_adjoint.end(),0.0);
    sys.v_adjoint[seed_dof] = seed;
    for (int i=0; i<steps; i++) {
      sys.update_state(dt,PassPhase::BackwardAdjoint);
    }
    grad_tau = global_field_by_name(sys,"dL_dparam_relaxation_time");
    grad_E = global_field_by_name(sys,"dL_dparam_youngs_modulus");
  }

  return { grad_tau, grad_E, phi };
}

struct ElementKickRun {
  Real grad_tau;
  Real grad_mu_e;
  Real phi;
};

ElementKickRun run_element_case(Real tau, Real mu_e, int steps, Real dt, Real seed, bool run_adjoint) {
  ElementSystem sys;

  const int Nnodes = 4;
  const int Ndofs_per_node = 2;
  const int Nelems = 1;
  const int Nnodes_per_elem = 4;

  double coordinates[8] = {
    0.0,0.0,
    1.0,0.0,
    1.0,1.0,
    0.0,1.0
  };
  double velocities[8] = {
    0.0,0.0,
    0.08,0.0,
    0.08,0.0,
    0.0,0.0
  };
  bool fixity[8] = {
    true,true,
    false,true,
    false,true,
    true,true
  };
  int connectivity[4] = { 0,1,2,3 };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 2.0;
  params["poissons_ratio"] = 0.25;
  params["relaxation_time"] = tau;
  params["shear_modulus_Maxwell_element"] = mu_e;
  params["mass_damping_factor"] = 0.0;
  params["contact_stiffness"] = 0.0;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity,Nelems,Nnodes_per_elem,params);
  sys.initialize_state();

  for (int i=0; i<steps; i++) {
    sys.update_state(dt,PassPhase::Forward);
  }

  const int seed_dof = 2; // node 1, x
  Real phi = seed*Real(sys.v[seed_dof].first);

  Real grad_tau = 0.0;
  Real grad_mu_e = 0.0;
  if (run_adjoint) {
    std::fill(sys.u_adjoint.begin(),sys.u_adjoint.end(),0.0);
    std::fill(sys.v_adjoint.begin(),sys.v_adjoint.end(),0.0);
    sys.v_adjoint[seed_dof] = seed;
    for (int i=0; i<steps; i++) {
      sys.update_state(dt,PassPhase::BackwardAdjoint);
    }
    grad_tau = global_field_by_name(sys,"dL_dparam_relaxation_time");
    grad_mu_e = global_field_by_name(sys,"dL_dparam_shear_modulus_Maxwell_element");
  }

  return { grad_tau, grad_mu_e, phi };
}

} // anonymous namespace

TEST(test_SystemAdjointKick, truss_delta_matches_terminal_velocity_fd) {
  const Real tau = 0.35;
  const Real E = 2.1;
  const int steps = 16;
  const Real dt = 2.0e-3;
  const Real seed = 0.7;

  TrussKickRun base = run_truss_case(tau,E,steps,dt,0.0,true);
  TrussKickRun seeded = run_truss_case(tau,E,steps,dt,seed,true);

  const Real h_tau = 1.0e-4;
  const Real h_E = 1.0e-4;

  Real phi_tau_plus = run_truss_case(tau+h_tau,E,steps,dt,seed,false).phi;
  Real phi_tau_minus = run_truss_case(tau-h_tau,E,steps,dt,seed,false).phi;
  Real fd_tau = (phi_tau_plus - phi_tau_minus)/(2.0*h_tau);

  Real phi_E_plus = run_truss_case(tau,E+h_E,steps,dt,seed,false).phi;
  Real phi_E_minus = run_truss_case(tau,E-h_E,steps,dt,seed,false).phi;
  Real fd_E = (phi_E_plus - phi_E_minus)/(2.0*h_E);

  Real delta_tau = seeded.grad_tau - base.grad_tau;
  Real delta_E = seeded.grad_E - base.grad_E;

  Real tol_tau = 8.0e-2 * std::max(Real(1.0),std::fabs(fd_tau));
  Real tol_E = 8.0e-2 * std::max(Real(1.0),std::fabs(fd_E));
  ASSERT_NEAR(delta_tau,fd_tau,tol_tau);
  ASSERT_NEAR(delta_E,fd_E,tol_E);
}

TEST(test_SystemAdjointKick, element_delta_matches_terminal_velocity_fd) {
  const Real tau = 0.33;
  const Real mu_e = 0.85;
  const int steps = 14;
  const Real dt = 1.5e-3;
  const Real seed = 0.5;

  ElementKickRun base = run_element_case(tau,mu_e,steps,dt,0.0,true);
  ElementKickRun seeded = run_element_case(tau,mu_e,steps,dt,seed,true);

  const Real h_tau = 8.0e-5;
  const Real h_mu = 8.0e-5;

  Real phi_tau_plus = run_element_case(tau+h_tau,mu_e,steps,dt,seed,false).phi;
  Real phi_tau_minus = run_element_case(tau-h_tau,mu_e,steps,dt,seed,false).phi;
  Real fd_tau = (phi_tau_plus - phi_tau_minus)/(2.0*h_tau);

  Real phi_mu_plus = run_element_case(tau,mu_e+h_mu,steps,dt,seed,false).phi;
  Real phi_mu_minus = run_element_case(tau,mu_e-h_mu,steps,dt,seed,false).phi;
  Real fd_mu = (phi_mu_plus - phi_mu_minus)/(2.0*h_mu);

  Real delta_tau = seeded.grad_tau - base.grad_tau;
  Real delta_mu = seeded.grad_mu_e - base.grad_mu_e;

  Real tol_tau = 1.2e-1 * std::max(Real(1.0),std::fabs(fd_tau));
  Real tol_mu = 1.2e-1 * std::max(Real(1.0),std::fabs(fd_mu));
  ASSERT_NEAR(delta_tau,fd_tau,tol_tau);
  ASSERT_NEAR(delta_mu,fd_mu,tol_mu);
}
