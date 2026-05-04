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
#include <array>

namespace {

using TrussSystem = System<
  Element<Material>,
  Truss<UniaxialViscoelasticity<Real,Real,Real> >,
  Real, Real, Real>;

using ElementSystem = System<
  Element<ViscoElasticity<Real,Real> >,
  Truss<UniaxialMaterial>,
  Real, Real, Real>;

using PointMassSystem = System<
  Element<Material>,
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

Real node_field_by_name(SystemBase& sys, int node_id, const std::string& name) {
  const int n_nodes = sys.get_num_entities("node");
  const int n_fields = sys.get_num_fields("node");
  if ((node_id < 0) || (node_id >= n_nodes)) { return 0.0; }

  std::vector<double> values(n_nodes*n_fields,0.0);
  sys.get_fields("node",values.data());
  for (int i=0; i<n_fields; i++) {
    if (name == sys.get_field_name("node",i)) {
      return values[n_fields*node_id + i];
    }
  }
  return 0.0;
}

Real element_field_by_name(SystemBase& sys, int elem_id, const std::string& name) {
  const int n_elem = sys.get_num_entities("element");
  const int n_fields = sys.get_num_fields("element");
  if ((elem_id < 0) || (elem_id >= n_elem)) { return 0.0; }

  std::vector<double> values(n_elem*n_fields,0.0);
  sys.get_fields("element",values.data());
  for (int i=0; i<n_fields; i++) {
    if (name == sys.get_field_name("element",i)) {
      return values[n_fields*elem_id + i];
    }
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

struct PointMassKickRun {
  Real grad_alpha;
  Real grad_k;
  Real phi;
};

PointMassKickRun run_point_mass_case(Real alpha_coeff, Real wall_k, Real y0, Real vy0,
                                     int steps, Real dt, Real seed, bool run_adjoint,
                                     bool use_seed_api = false) {
  PointMassSystem sys;

  const int Nnodes = 1;
  const int Ndofs_per_node = 2;
  const int Nelems = 0;
  const int Nnodes_per_elem = 4;

  double coordinates[2] = { 0.0, y0 };
  double velocities[2]  = { 0.0, vy0 };
  bool fixity[2]        = { true, false };
  int connectivity_dummy[4] = { 0,0,0,0 };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;
  params["mass_damping_factor"] = alpha_coeff;
  params["contact_stiffness"] = wall_k;
  params["body_force_x"] = 0.0;
  params["body_force_y"] = 0.0;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity_dummy,Nelems,Nnodes_per_elem,params);
  int point_ids[1] = { 0 };
  double point_mass[1] = { 1.0 };
  sys.initialize_point_mass(point_ids,point_mass,1,params);
  sys.initialize_state();

  for (int i=0; i<steps; i++) {
    sys.update_state(dt,PassPhase::Forward);
  }

  const int seed_dof = 1; // y dof
  Real phi = seed*Real(sys.v[seed_dof].first);

  Real grad_alpha = 0.0;
  Real grad_k = 0.0;
  if (run_adjoint) {
    if (use_seed_api) {
      sys.clear_adjoint_state();
      const int sensor_nodes[1] = { 0 };
      const double seed_xy[2] = { 0.0, seed };
      sys.add_nodal_velocity_adjoint_seed(sensor_nodes,seed_xy,1);
      const double zero_disp_seed[2] = { 0.0, 0.0 };
      sys.add_nodal_displacement_adjoint_seed(sensor_nodes,zero_disp_seed,1);
    } else {
      std::fill(sys.u_adjoint.begin(),sys.u_adjoint.end(),0.0);
      std::fill(sys.v_adjoint.begin(),sys.v_adjoint.end(),0.0);
      sys.v_adjoint[seed_dof] = seed;
    }
    for (int i=0; i<steps; i++) {
      sys.update_state(dt,PassPhase::BackwardAdjoint);
    }
    grad_alpha = global_field_by_name(sys,"dL_dparam_mass_damping_factor");
    grad_k = global_field_by_name(sys,"dL_dparam_contact_stiffness");
  }

  return { grad_alpha, grad_k, phi };
}

Real two_layer_scale_coeffs[2] = { 1.0, 1.0 };
Real two_layer_tau_coeffs[2] = { 0.1, 0.1 };

double two_layer_stiffness_scaling(double, double y) {
  return (y < 1.0) ? two_layer_scale_coeffs[0] : two_layer_scale_coeffs[1];
}

double two_layer_relaxation_time(double, double y) {
  return (y < 1.0) ? two_layer_tau_coeffs[0] : two_layer_tau_coeffs[1];
}

struct TwoLayerInverseRun {
  Real loss;
  Real grad_tau;
  std::array<Real,2> grad_layers;
  std::vector<Real> sensor_history; // [vy_sensor0, vy_sensor1] per time step
};

TwoLayerInverseRun run_two_layer_inverse_case(Real scale0, Real scale1, Real tau,
                                              int steps, Real dt,
                                              const std::vector<Real>* observed_history,
                                              bool compute_gradients) {
  ElementSystem sys;

  const int Nnodes = 6;
  const int Ndofs_per_node = 2;
  const int Nelems = 2;
  const int Nnodes_per_elem = 4;

  double coordinates[12] = {
    0.0,0.0, 1.0,0.0,
    0.0,1.0, 1.0,1.0,
    0.0,2.0, 1.0,2.0
  };
  double velocities[12] = {
    0.0,0.0, 0.0,0.0,
    0.0,0.0, 0.0,0.0,
    0.0,-0.18, 0.0,-0.18
  };
  bool fixity[12] = {
    true,true,  true,true,   // bottom row fixed
    true,false, true,false,  // middle row: x fixed, y free
    true,false, true,false   // top row: x fixed, y free
  };
  int connectivity[8] = {
    0,1,3,2,
    2,3,5,4
  };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 2.0;
  params["poissons_ratio"] = 0.25;
  params["relaxation_time"] = tau;
  params["shear_modulus_Maxwell_element"] = 0.8;
  params["mass_damping_factor"] = 0.0;
  params["contact_stiffness"] = 0.0;
  params["adjoint_material_objective_weight"] = 0.0; // sensor-misfit-only objective

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity,Nelems,Nnodes_per_elem,params);
  two_layer_scale_coeffs[0] = scale0;
  two_layer_scale_coeffs[1] = scale1;
  sys.initialize_variable_properties(two_layer_stiffness_scaling);
  sys.initialize_state();

  std::vector<Real> sensor_history(2*steps,0.0);
  Real loss = 0.0;
  for (int k=0; k<steps; k++) {
    sys.update_state(dt,PassPhase::Forward);
    const Real vy0 = Real(sys.v[2*4 + 1].first);
    const Real vy1 = Real(sys.v[2*5 + 1].first);
    sensor_history[2*k + 0] = vy0;
    sensor_history[2*k + 1] = vy1;
    if (observed_history != nullptr) {
      const Real r0 = vy0 - (*observed_history)[2*k + 0];
      const Real r1 = vy1 - (*observed_history)[2*k + 1];
      loss += 0.5*(r0*r0 + r1*r1);
    }
  }

  Real grad_tau = 0.0;
  std::array<Real,2> grad_layers = { 0.0, 0.0 };
  if (compute_gradients && (observed_history != nullptr)) {
    sys.clear_adjoint_state();
    const int sensor_nodes[2] = { 4, 5 };
    for (int rev=0; rev<steps; rev++) {
      const int k = steps - 1 - rev;
      const double seeds_xy[4] = {
        0.0, sensor_history[2*k + 0] - (*observed_history)[2*k + 0],
        0.0, sensor_history[2*k + 1] - (*observed_history)[2*k + 1]
      };
      sys.add_nodal_velocity_adjoint_seed(sensor_nodes,seeds_xy,2);
      sys.update_state(dt,PassPhase::BackwardAdjoint);
    }

    grad_tau = global_field_by_name(sys,"dL_dparam_relaxation_time");

    int stiffness_grad_field = -1;
    const int num_elem_fields = sys.get_num_fields("element");
    for (int i=0; i<num_elem_fields; i++) {
      const std::string field_name = sys.get_field_name("element",i);
      if (field_name == "dparam_stiffness_scaling_factor") {
        stiffness_grad_field = i;
        break;
      }
    }
    if (stiffness_grad_field < 0) {
      std::cout << "Missing element field: dparam_stiffness_scaling_factor" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::vector<double> element_fields(Nelems*num_elem_fields,0.0);
    sys.get_fields("element",element_fields.data());
    grad_layers[0] = element_fields[num_elem_fields*0 + stiffness_grad_field];
    grad_layers[1] = element_fields[num_elem_fields*1 + stiffness_grad_field];
  }

  return { loss, grad_tau, grad_layers, sensor_history };
}

struct TwoLayerSpatialTauInverseRun {
  Real loss;
  Real grad_tau_sum;
  std::array<Real,2> grad_tau_by_element;
  std::array<Real,2> tau_by_element;
  std::vector<Real> sensor_history; // [vy_sensor0, vy_sensor1] per time step
};

TwoLayerSpatialTauInverseRun run_two_layer_spatial_tau_inverse_case(
    Real scale0, Real scale1, Real tau0, Real tau1,
    int steps, Real dt,
    const std::vector<Real>* observed_history,
    bool compute_gradients) {
  ElementSystem sys;

  const int Nnodes = 6;
  const int Ndofs_per_node = 2;
  const int Nelems = 2;
  const int Nnodes_per_elem = 4;

  double coordinates[12] = {
    0.0,0.0, 1.0,0.0,
    0.0,1.0, 1.0,1.0,
    0.0,2.0, 1.0,2.0
  };
  double velocities[12] = {
    0.0,0.0, 0.0,0.0,
    0.0,0.0, 0.0,0.0,
    0.0,-0.18, 0.0,-0.18
  };
  bool fixity[12] = {
    true,true,  true,true,
    true,false, true,false,
    true,false, true,false
  };
  int connectivity[8] = {
    0,1,3,2,
    2,3,5,4
  };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 2.0;
  params["poissons_ratio"] = 0.25;
  params["relaxation_time"] = 0.1; // fallback scalar; overridden by spatial callback
  params["shear_modulus_Maxwell_element"] = 0.8;
  params["mass_damping_factor"] = 0.0;
  params["contact_stiffness"] = 0.0;
  params["adjoint_material_objective_weight"] = 0.0;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity,Nelems,Nnodes_per_elem,params);

  two_layer_scale_coeffs[0] = scale0;
  two_layer_scale_coeffs[1] = scale1;
  two_layer_tau_coeffs[0] = tau0;
  two_layer_tau_coeffs[1] = tau1;
  sys.initialize_variable_properties(two_layer_stiffness_scaling);
  sys.initialize_variable_relaxation_time(two_layer_relaxation_time);
  sys.initialize_state();

  std::vector<Real> sensor_history(2*steps,0.0);
  Real loss = 0.0;
  for (int k=0; k<steps; k++) {
    sys.update_state(dt,PassPhase::Forward);
    const Real vy0 = Real(sys.v[2*4 + 1].first);
    const Real vy1 = Real(sys.v[2*5 + 1].first);
    sensor_history[2*k + 0] = vy0;
    sensor_history[2*k + 1] = vy1;
    if (observed_history != nullptr) {
      const Real r0 = vy0 - (*observed_history)[2*k + 0];
      const Real r1 = vy1 - (*observed_history)[2*k + 1];
      loss += 0.5*(r0*r0 + r1*r1);
    }
  }

  const Real tau_e0 = element_field_by_name(sys,0,"relaxation_time_local");
  const Real tau_e1 = element_field_by_name(sys,1,"relaxation_time_local");

  Real grad_tau_sum = 0.0;
  std::array<Real,2> grad_tau_by_element = { 0.0, 0.0 };
  if (compute_gradients && (observed_history != nullptr)) {
    sys.clear_adjoint_state();
    const int sensor_nodes[2] = { 4, 5 };
    for (int rev=0; rev<steps; rev++) {
      const int k = steps - 1 - rev;
      const double seeds_xy[4] = {
        0.0, sensor_history[2*k + 0] - (*observed_history)[2*k + 0],
        0.0, sensor_history[2*k + 1] - (*observed_history)[2*k + 1]
      };
      sys.add_nodal_velocity_adjoint_seed(sensor_nodes,seeds_xy,2);
      sys.update_state(dt,PassPhase::BackwardAdjoint);
    }

    grad_tau_sum = global_field_by_name(sys,"dL_dparam_relaxation_time");
    grad_tau_by_element[0] = element_field_by_name(sys,0,"dparam_relaxation_time");
    grad_tau_by_element[1] = element_field_by_name(sys,1,"dparam_relaxation_time");
  }

  return { loss, grad_tau_sum, grad_tau_by_element, { tau_e0, tau_e1 }, sensor_history };
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

TEST(test_SystemAdjointKick, damping_only_delta_matches_terminal_velocity_fd) {
  const Real alpha_coeff = 0.35;
  const Real wall_k = 0.0;
  const Real y0 = 0.25;
  const Real vy0 = 0.18;
  const int steps = 30;
  const Real dt = 2.0e-3;
  const Real seed = 0.8;

  PointMassKickRun base = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,0.0,true);
  PointMassKickRun seeded = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,seed,true);

  const Real h_alpha = 1.0e-4;
  Real phi_plus = run_point_mass_case(alpha_coeff+h_alpha,wall_k,y0,vy0,steps,dt,seed,false).phi;
  Real phi_minus = run_point_mass_case(alpha_coeff-h_alpha,wall_k,y0,vy0,steps,dt,seed,false).phi;
  Real fd_alpha = (phi_plus - phi_minus)/(2.0*h_alpha);

  Real delta_alpha = seeded.grad_alpha - base.grad_alpha;
  Real delta_k = seeded.grad_k - base.grad_k;

  Real tol_alpha = 8.0e-2 * std::max(Real(1.0),std::fabs(fd_alpha));
  ASSERT_NEAR(delta_alpha,fd_alpha,tol_alpha);
  ASSERT_NEAR(delta_k,0.0,1.0e-12);
}

TEST(test_SystemAdjointKick, rigid_wall_only_delta_matches_terminal_velocity_fd) {
  const Real alpha_coeff = 0.0;
  const Real wall_k = 7.0;
  const Real y0 = -0.3;
  const Real vy0 = 0.0;
  const int steps = 24;
  const Real dt = 1.0e-3;
  const Real seed = 0.6;

  PointMassKickRun base = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,0.0,true);
  PointMassKickRun seeded = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,seed,true);

  const Real h_k = 1.0e-4;
  Real phi_plus = run_point_mass_case(alpha_coeff,wall_k+h_k,y0,vy0,steps,dt,seed,false).phi;
  Real phi_minus = run_point_mass_case(alpha_coeff,wall_k-h_k,y0,vy0,steps,dt,seed,false).phi;
  Real fd_k = (phi_plus - phi_minus)/(2.0*h_k);

  Real delta_alpha = seeded.grad_alpha - base.grad_alpha;
  Real delta_k = seeded.grad_k - base.grad_k;

  Real tol_k = 1.2e-1 * std::max(Real(1.0),std::fabs(fd_k));
  ASSERT_NEAR(delta_k,fd_k,tol_k);
  ASSERT_NEAR(delta_alpha,0.0,1.0e-12);
}

TEST(test_SystemAdjointKick, damping_and_wall_delta_match_terminal_velocity_fd) {
  const Real alpha_coeff = 0.22;
  const Real wall_k = 5.0;
  const Real y0 = -0.3;
  const Real vy0 = 0.05;
  const int steps = 24;
  const Real dt = 1.0e-3;
  const Real seed = 0.65;

  PointMassKickRun base = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,0.0,true);
  PointMassKickRun seeded = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,seed,true);

  const Real h_alpha = 1.0e-4;
  const Real h_k = 1.0e-4;

  Real phi_alpha_plus = run_point_mass_case(alpha_coeff+h_alpha,wall_k,y0,vy0,steps,dt,seed,false).phi;
  Real phi_alpha_minus = run_point_mass_case(alpha_coeff-h_alpha,wall_k,y0,vy0,steps,dt,seed,false).phi;
  Real fd_alpha = (phi_alpha_plus - phi_alpha_minus)/(2.0*h_alpha);

  Real phi_k_plus = run_point_mass_case(alpha_coeff,wall_k+h_k,y0,vy0,steps,dt,seed,false).phi;
  Real phi_k_minus = run_point_mass_case(alpha_coeff,wall_k-h_k,y0,vy0,steps,dt,seed,false).phi;
  Real fd_k = (phi_k_plus - phi_k_minus)/(2.0*h_k);

  Real delta_alpha = seeded.grad_alpha - base.grad_alpha;
  Real delta_k = seeded.grad_k - base.grad_k;

  Real tol_alpha = 1.5e-1 * std::max(Real(1.0),std::fabs(fd_alpha));
  Real tol_k = 1.5e-1 * std::max(Real(1.0),std::fabs(fd_k));
  ASSERT_NEAR(delta_alpha,fd_alpha,tol_alpha);
  ASSERT_NEAR(delta_k,fd_k,tol_k);
}

TEST(test_SystemAdjointKick, nodal_velocity_seed_api_matches_direct_seed) {
  const Real alpha_coeff = 0.17;
  const Real wall_k = 2.8;
  const Real y0 = -0.2;
  const Real vy0 = 0.11;
  const int steps = 20;
  const Real dt = 1.25e-3;
  const Real seed = 0.55;

  PointMassKickRun direct = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,seed,true,false);
  PointMassKickRun api = run_point_mass_case(alpha_coeff,wall_k,y0,vy0,steps,dt,seed,true,true);

  ASSERT_NEAR(api.phi,direct.phi,1.0e-14);
  ASSERT_NEAR(api.grad_alpha,direct.grad_alpha,1.0e-12);
  ASSERT_NEAR(api.grad_k,direct.grad_k,1.0e-12);
}

TEST(test_SystemAdjointKick, node_fields_expose_global_adjoint_arrays) {
  PointMassSystem sys;

  const int Nnodes = 1;
  const int Ndofs_per_node = 2;
  const int Nelems = 0;
  const int Nnodes_per_elem = 4;

  double coordinates[2] = { 0.0, 0.0 };
  double velocities[2]  = { 0.0, 0.0 };
  bool fixity[2]        = { false, false };
  int connectivity_dummy[4] = { 0,0,0,0 };
  int point_ids[1] = { 0 };
  double point_mass[1] = { 1.0 };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity_dummy,Nelems,Nnodes_per_elem,params);
  sys.initialize_point_mass(point_ids,point_mass,1,params);
  sys.initialize_state();

  sys.clear_adjoint_state();
  const int sensor_nodes[1] = { 0 };
  const double velocity_seed_xy[2] = { 0.3, -0.7 };
  const double displacement_seed_xy[2] = { -0.2, 0.5 };
  sys.add_nodal_velocity_adjoint_seed(sensor_nodes,velocity_seed_xy,1);
  sys.add_nodal_displacement_adjoint_seed(sensor_nodes,displacement_seed_xy,1);

  ASSERT_NEAR(node_field_by_name(sys,0,"adjoint_velocity_X"),0.3,1.0e-14);
  ASSERT_NEAR(node_field_by_name(sys,0,"adjoint_velocity_Y"),-0.7,1.0e-14);
  ASSERT_NEAR(node_field_by_name(sys,0,"adjoint_displacement_X"),-0.2,1.0e-14);
  ASSERT_NEAR(node_field_by_name(sys,0,"adjoint_displacement_Y"),0.5,1.0e-14);
}

TEST(test_SystemAdjointKick, invalid_sensor_node_id_is_rejected) {
  PointMassSystem sys;

  const int Nnodes = 1;
  const int Ndofs_per_node = 2;
  const int Nelems = 0;
  const int Nnodes_per_elem = 4;

  double coordinates[2] = { 0.0, 0.0 };
  double velocities[2]  = { 0.0, 0.0 };
  bool fixity[2]        = { false, false };
  int connectivity_dummy[4] = { 0,0,0,0 };
  int point_ids[1] = { 0 };
  double point_mass[1] = { 1.0 };

  Parameters params;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;

  sys.initialize(coordinates,velocities,fixity,Nnodes,Ndofs_per_node,
                 connectivity_dummy,Nelems,Nnodes_per_elem,params);
  sys.initialize_point_mass(point_ids,point_mass,1,params);
  sys.initialize_state();

  const int bad_nodes[1] = { 5 };
  const double seeds[2] = { 0.0, 1.0 };
  ASSERT_DEATH(sys.add_nodal_velocity_adjoint_seed(bad_nodes,seeds,1), "");
}

TEST(test_SystemAdjointKick, two_layer_sensor_misfit_gradients_match_finite_difference) {
  const int steps = 14;
  const Real dt = 1.0e-3;

  const Real true_scale0 = 0.78;
  const Real true_scale1 = 1.22;
  const Real true_tau = 0.14;
  TwoLayerInverseRun observed = run_two_layer_inverse_case(
    true_scale0,true_scale1,true_tau,steps,dt,nullptr,false
  );

  const Real cand_scale0 = 1.05;
  const Real cand_scale1 = 0.92;
  const Real cand_tau = 0.22;
  TwoLayerInverseRun adj = run_two_layer_inverse_case(
    cand_scale0,cand_scale1,cand_tau,steps,dt,&observed.sensor_history,true
  );

  const Real h_tau = 1.0e-4;
  const Real h_scale = 1.0e-4;

  const Real loss_tau_plus = run_two_layer_inverse_case(
    cand_scale0,cand_scale1,cand_tau + h_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real loss_tau_minus = run_two_layer_inverse_case(
    cand_scale0,cand_scale1,cand_tau - h_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real fd_tau = (loss_tau_plus - loss_tau_minus)/(2.0*h_tau);

  const Real loss_s0_plus = run_two_layer_inverse_case(
    cand_scale0 + h_scale,cand_scale1,cand_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real loss_s0_minus = run_two_layer_inverse_case(
    cand_scale0 - h_scale,cand_scale1,cand_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real fd_s0 = (loss_s0_plus - loss_s0_minus)/(2.0*h_scale);

  const Real loss_s1_plus = run_two_layer_inverse_case(
    cand_scale0,cand_scale1 + h_scale,cand_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real loss_s1_minus = run_two_layer_inverse_case(
    cand_scale0,cand_scale1 - h_scale,cand_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real fd_s1 = (loss_s1_plus - loss_s1_minus)/(2.0*h_scale);

  auto rel_err = [](Real adj_grad, Real fd_grad) {
    return std::fabs(adj_grad - fd_grad) /
      std::max({std::fabs(adj_grad),std::fabs(fd_grad),Real(1.0e-12)});
  };

  ASSERT_LT(rel_err(adj.grad_tau,fd_tau),1.0e-1);
  ASSERT_LT(rel_err(adj.grad_layers[0],fd_s0),1.0e-1);
  ASSERT_LT(rel_err(adj.grad_layers[1],fd_s1),1.0e-1);
}

TEST(test_SystemAdjointKick, spatial_tau_initialization_sets_element_fields) {
  const int steps = 4;
  const Real dt = 1.0e-3;

  const Real scale0 = 0.95;
  const Real scale1 = 1.05;
  const Real tau0 = 0.07;
  const Real tau1 = 0.16;

  TwoLayerSpatialTauInverseRun run = run_two_layer_spatial_tau_inverse_case(
    scale0,scale1,tau0,tau1,steps,dt,nullptr,false
  );

  ASSERT_NEAR(run.tau_by_element[0],tau0,1.0e-12);
  ASSERT_NEAR(run.tau_by_element[1],tau1,1.0e-12);
}

TEST(test_SystemAdjointKick, spatial_tau_element_gradients_match_finite_difference) {
  const int steps = 14;
  const Real dt = 1.0e-3;

  const Real true_scale0 = 0.82;
  const Real true_scale1 = 1.18;
  const Real true_tau0 = 0.08;
  const Real true_tau1 = 0.15;
  TwoLayerSpatialTauInverseRun observed = run_two_layer_spatial_tau_inverse_case(
    true_scale0,true_scale1,true_tau0,true_tau1,steps,dt,nullptr,false
  );

  const Real cand_scale0 = 1.06;
  const Real cand_scale1 = 0.91;
  const Real cand_tau0 = 0.20;
  const Real cand_tau1 = 0.11;
  TwoLayerSpatialTauInverseRun adj = run_two_layer_spatial_tau_inverse_case(
    cand_scale0,cand_scale1,cand_tau0,cand_tau1,steps,dt,&observed.sensor_history,true
  );

  const Real h_tau = 1.0e-4;
  const Real loss_t0_plus = run_two_layer_spatial_tau_inverse_case(
    cand_scale0,cand_scale1,cand_tau0 + h_tau,cand_tau1,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real loss_t0_minus = run_two_layer_spatial_tau_inverse_case(
    cand_scale0,cand_scale1,cand_tau0 - h_tau,cand_tau1,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real fd_t0 = (loss_t0_plus - loss_t0_minus)/(2.0*h_tau);

  const Real loss_t1_plus = run_two_layer_spatial_tau_inverse_case(
    cand_scale0,cand_scale1,cand_tau0,cand_tau1 + h_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real loss_t1_minus = run_two_layer_spatial_tau_inverse_case(
    cand_scale0,cand_scale1,cand_tau0,cand_tau1 - h_tau,steps,dt,&observed.sensor_history,false
  ).loss;
  const Real fd_t1 = (loss_t1_plus - loss_t1_minus)/(2.0*h_tau);

  auto rel_err = [](Real adj_grad, Real fd_grad) {
    return std::fabs(adj_grad - fd_grad) /
      std::max({std::fabs(adj_grad),std::fabs(fd_grad),Real(1.0e-12)});
  };

  ASSERT_LT(rel_err(adj.grad_tau_by_element[0],fd_t0),1.2e-1);
  ASSERT_LT(rel_err(adj.grad_tau_by_element[1],fd_t1),1.2e-1);
  ASSERT_NEAR(adj.grad_tau_sum,adj.grad_tau_by_element[0] + adj.grad_tau_by_element[1],1.0e-10);
}
