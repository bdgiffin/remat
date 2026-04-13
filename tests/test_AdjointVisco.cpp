#include "gtest/gtest.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "Truss.h"
#include "UniaxialViscoelasticity.h"
#include "Fixed.h"
#include "Rational.h"
#include "Parameters.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <string>

namespace {

using FloatTrussVisco = Truss<UniaxialViscoelasticity<Real,Real,Real>>;
using FixedAdjFloatTrussVisco = Truss<UniaxialViscoelasticity<Fixed_E,Rational,Real>>;

using AdjointFloatSystem = System<Element<Material>,FloatTrussVisco,Real,Real,Real,true>;
using AdjointFixedAdjFloatSystem = System<Element<Material>,FixedAdjFloatTrussVisco,Fixed_V,Fixed_U,Rational,true>;

struct RunResult {
  Real loss = 0.0;
  Real dL_drelaxation_time_sum_truss = 0.0;
  Real dL_drelaxation_time_global = 0.0;
};

template<class SystemT>
Real get_global_field(SystemT& problem, const std::string& field_name) {
  const int Nglobal = problem.get_num_fields("global");
  std::vector<double> global_fields(Nglobal,0.0);
  problem.get_fields("global",global_fields.data());
  for (int i=0; i<Nglobal; ++i) {
    const char* current_name = problem.get_field_name("global",i);
    if (current_name != nullptr && field_name == current_name) {
      return Real(global_fields[i]);
    }
  }
  return 0.0;
}

Parameters make_params(Real tau, Real E, int mat_overflow_limit) {
  Parameters params;
  params["body_force_x"] = 0.0;
  params["body_force_y"] = 0.0;
  params["mass_damping_factor"] = 0.0;
  params["dt_scale_factor"] = 1.0;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;
  params["truss_density"] = 1.0;
  params["truss_youngs_modulus"] = E;
  params["area"] = 1.0;
  params["relaxation_time"] = tau;
  params["mat_overflow_limit"] = mat_overflow_limit;
  params["overflow_limit"] = 1000000.0;
  return params;
}

template<class SystemT>
RunResult run_chain_case(Real tau, Real E, int num_elements, Real dt, int nsteps, Real impact_velocity, bool run_backward) {
  const int Nnodes = num_elements + 1;
  const int Ndofs_per_node = 2;
  const int Ndofs = Nnodes*Ndofs_per_node;

  std::vector<double> coordinates(Ndofs,0.0);
  std::vector<double> velocities(Ndofs,0.0);
  std::unique_ptr<bool[]> fixity(new bool[Ndofs]);
  std::vector<int> connectivity(1,0); // no solid elements in this test
  std::vector<int> truss_connectivity(2*num_elements,0);

  for (int i=0; i<Ndofs; ++i) fixity[i] = false;

  for (int i=0; i<Nnodes; ++i) {
    coordinates[2*i+0] = Real(i);
    coordinates[2*i+1] = 0.0;
    fixity[2*i+1] = true;
  }
  fixity[0] = true;
  velocities[2*(Nnodes-1)+0] = impact_velocity;

  for (int e=0; e<num_elements; ++e) {
    truss_connectivity[2*e+0] = e;
    truss_connectivity[2*e+1] = e+1;
  }

  Parameters params = make_params(tau,E,1000000);

  SystemT problem;
  problem.initialize(coordinates.data(),velocities.data(),fixity.get(),Nnodes,Ndofs_per_node,
                     connectivity.data(),0,4,params);
  problem.initialize_truss_elements(truss_connectivity.data(),num_elements,params);
  problem.initialize_state();

  RunResult result;
  const int Nstate_vars_per_truss = problem.m_truss.num_state_vars();

  for (int step=0; step<nsteps; ++step) {
    for (int e=0; e<num_elements; ++e) {
      const Real sigma = problem.truss_state[Nstate_vars_per_truss*e+0];
      result.loss += 0.5*sigma*sigma/E;
    }
    problem.update_state(+dt);
  }

  if (run_backward) {
    for (int step=0; step<nsteps; ++step) {
      problem.update_state(-dt);
    }
  }

  for (int e=0; e<num_elements; ++e) {
    result.dL_drelaxation_time_sum_truss += problem.truss_state[Nstate_vars_per_truss*e+5];
  }
  result.dL_drelaxation_time_global = get_global_field(problem,"dL_dparam_relaxation_time");

  return result;
}

template<class SystemT>
void check_primal_reversibility(int num_elements, Real tau, Real E, Real dt, int nsteps, Real impact_velocity, Real tol) {
  const int Nnodes = num_elements + 1;
  const int Ndofs_per_node = 2;
  const int Ndofs = Nnodes*Ndofs_per_node;

  std::vector<double> coordinates(Ndofs,0.0);
  std::vector<double> velocities(Ndofs,0.0);
  std::unique_ptr<bool[]> fixity(new bool[Ndofs]);
  std::vector<int> connectivity(1,0); // no solid elements in this test
  std::vector<int> truss_connectivity(2*num_elements,0);

  for (int i=0; i<Ndofs; ++i) fixity[i] = false;

  for (int i=0; i<Nnodes; ++i) {
    coordinates[2*i+0] = Real(i);
    coordinates[2*i+1] = 0.0;
    fixity[2*i+1] = true;
  }
  fixity[0] = true;
  velocities[2*(Nnodes-1)+0] = impact_velocity;

  for (int e=0; e<num_elements; ++e) {
    truss_connectivity[2*e+0] = e;
    truss_connectivity[2*e+1] = e+1;
  }

  Parameters params = make_params(tau,E,1000000);
  SystemT problem;
  problem.initialize(coordinates.data(),velocities.data(),fixity.get(),Nnodes,Ndofs_per_node,
                     connectivity.data(),0,4,params);
  problem.initialize_truss_elements(truss_connectivity.data(),num_elements,params);
  problem.initialize_state();

  std::vector<Real> u_initial(Ndofs,0.0);
  std::vector<Real> v_initial(Ndofs,0.0);
  for (int i=0; i<Ndofs; ++i) {
    u_initial[i] = Real(problem.u[i].first);
    v_initial[i] = Real(problem.v[i].first);
  }

  const int Nstate_vars_per_truss = problem.m_truss.num_state_vars();
  std::vector<Real> truss_initial(5*num_elements,0.0);
  for (int e=0; e<num_elements; ++e) {
    for (int s=0; s<5; ++s) {
      truss_initial[5*e+s] = problem.truss_state[Nstate_vars_per_truss*e+s];
    }
  }

  for (int step=0; step<nsteps; ++step) {
    problem.update_state(+dt);
  }
  for (int step=0; step<nsteps; ++step) {
    problem.update_state(-dt);
  }

  for (int i=0; i<Ndofs; ++i) {
    ASSERT_NEAR(Real(problem.u[i].first),u_initial[i],tol);
    ASSERT_NEAR(Real(problem.v[i].first),v_initial[i],tol);
  }

  for (int e=0; e<num_elements; ++e) {
    for (int s=0; s<5; ++s) {
      ASSERT_NEAR(problem.truss_state[Nstate_vars_per_truss*e+s],truss_initial[5*e+s],tol);
    }
  }
}

template<class SystemT>
void check_gradient_matches_fd(int num_elements, Real dt, int nsteps, Real tau, Real E, Real impact_velocity, Real h_tau, Real rel_tol) {
  RunResult adj = run_chain_case<SystemT>(tau,E,num_elements,dt,nsteps,impact_velocity,true);
  RunResult plus = run_chain_case<SystemT>(tau+h_tau,E,num_elements,dt,nsteps,impact_velocity,false);
  RunResult minus = run_chain_case<SystemT>(tau-h_tau,E,num_elements,dt,nsteps,impact_velocity,false);

  const Real fd_tau = (plus.loss - minus.loss)/(2.0*h_tau);
  const Real tol = rel_tol*std::max(Real(1.0),std::fabs(fd_tau));
  ASSERT_NEAR(adj.dL_drelaxation_time_global,adj.dL_drelaxation_time_sum_truss,1.0e-12);
  ASSERT_NEAR(adj.dL_drelaxation_time_global,fd_tau,tol);
}

} // namespace

TEST(test_AdjointVisco, gradient_one_element_float_mode) {
  check_gradient_matches_fd<AdjointFloatSystem>(1,2.0e-3,220,0.35,2.0,0.25,1.0e-4,3.0e-2);
}

TEST(test_AdjointVisco, gradient_multi_element_float_mode) {
  check_gradient_matches_fd<AdjointFloatSystem>(4,1.0e-3,260,0.30,2.0,0.20,1.0e-4,5.0e-2);
}

TEST(test_AdjointVisco, gradient_one_element_fixed_adj_float_mode) {
  check_gradient_matches_fd<AdjointFixedAdjFloatSystem>(1,2.0e-3,220,0.35,2.0,0.25,2.0e-3,8.0e-2);
}

TEST(test_AdjointVisco, gradient_multi_element_fixed_adj_float_mode) {
  check_gradient_matches_fd<AdjointFixedAdjFloatSystem>(4,1.0e-3,260,0.30,2.0,0.20,2.0e-3,1.5e-1);
}

TEST(test_AdjointVisco, reversibility_float_mode) {
  check_primal_reversibility<AdjointFloatSystem>(3,0.25,2.0,1.0e-3,200,0.20,3.0e-4);
}

TEST(test_AdjointVisco, reversibility_fixed_adj_float_mode) {
  check_primal_reversibility<AdjointFixedAdjFloatSystem>(3,0.25,2.0,1.0e-3,200,0.20,6.0e-4);
}
