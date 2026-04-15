#ifndef TRUSS_H
#define TRUSS_H

#include "PassPhase.h"
#include "MaterialAdjoint.h"
#include "Parameters.h"
#include <vector>
#include <string>
#include <math.h>
#include <type_traits>
#include <utility>

template<typename Material_T>
class Truss {
private:
  template<typename M = Material_T>
  auto call_material_update(Real lambda, Real &psi, Real* state, Real dt, PassPhase phase, int)
    -> decltype(std::declval<M&>().update(lambda,psi,state,dt,phase), void()) {
    m_model.update(lambda,psi,state,dt,phase);
  }

  template<typename M = Material_T>
  void call_material_update(Real lambda, Real &psi, Real* state, Real dt, PassPhase phase, long) {
    Real signed_dt = is_reverse_phase(phase) ? -std::fabs(dt) : std::fabs(dt);
    m_model.update(lambda,psi,state,signed_dt);
  }

public:
  Material_T m_model; // material model

  // Empty constructor
  Truss(void) { }

  // Parameterized constructor
  Truss(Parameters& params) : m_model(params) { }
    
  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return m_model.num_state_vars(); }

  // Initialize the element state
  void initialize(Real* state) {
    
    // zero-initialize element-averaged state
    m_model.initialize(state);
    
  } // initialize()
    
  // Update the element state using the current nodal displacements
  void update(Real (&x)[4], Real (&u)[4], Real (&m)[4], Real (&f)[4], Real &E, Real* state, Real dt, PassPhase phase) {

    // Determine nodal masses:
    {
      // Compute mid-point shape function partial derivatives (dphi/dxi = 0.5*xi)
      const Real dxi[2] = { -0.5, +0.5 };

      // Compute mid-point Jacobian (J = dx/dxi = sum_a x_a (x) dphi/dxi)
      Real J[2] = { 0.0, 0.0 };
      for (int i = 0; i < num_nodes(); i++) {
	J[0] += x[2*i+0] * dxi[i];
	J[1] += x[2*i+1] * dxi[i];
      }

      // Compute norm of J
      Real normJ = std::sqrt(J[0]*J[0] + J[1]*J[1]);

      // Compute total element mass going to each node = (2 x normJ x mass_per_unit_length) / 2
      Real mass = normJ*m_model.mass_per_unit_length();

      // Sum nodal mass contributions
      for (int i = 0; i < 2*num_nodes(); i++) m[i] = mass;
    }

    // Determine nodal forces:
    {
      
      // Define parent nodal coordinates
      const Real xi[2] = { -1.0, +1.0 };

      // Set quadrature point position and weight
      const Real xiq = 0.0;
      const Real wq  = 2.0;
	
      // Compute shape function partial derivatives (dphi/dxi)
      const Real dxi[2] = { 0.5f*xi[0], 0.5f*xi[1] };

      // Compute initial and current Jacobians (J = dx/dxi = sum_a x_a (x) dphi/dxi)
      Real J0[2] = {0.0,0.0};
      Real J[2]  = {0.0,0.0};
      for (int i = 0; i < 2; i++) {
	J0[0] += x[2*i+0] * dxi[i];
	J0[1] += x[2*i+1] * dxi[i];
	Real xt = x[2*i+0] + u[2*i+0]; Real yt = x[2*i+1] + u[2*i+1];
	J[0] += xt * dxi[i];
	J[1] += yt * dxi[i];
      }

      // Compute norm of J0 and J
      Real normJ0 = std::sqrt(J0[0]*J0[0] + J0[1]*J0[1]);
      Real normJ  = std::sqrt(J[0]*J[0] + J[1]*J[1]);
	
      // Compute the current axial stretch ratio lambda = normJ/normJ0
      Real lambda = normJ/normJ0;

      // update the material state
      Real psi;
      call_material_update(lambda,psi,state,dt,phase,0);

      // Compute the tangent vector
      Real t[2] = { J[0]/normJ, J[1]/normJ };

      // Compute the axial force vector (state[0] = dpsi/dlambda)
      Real f_hat[2] = { state[0]*t[0], state[0]*t[1] };

      // sum nodal forces by stress divergence: f_hat.dphi/dxi * d[]
      for (int i = 0; i < num_nodes(); i++) {
	f[2*i+0] += f_hat[0]*dxi[i]*wq;
	f[2*i+1] += f_hat[1]*dxi[i]*wq;
      }

      // sum contribution to the total elastic strain energy
      E += psi*normJ0*wq;
	
    }

  } // update()

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state)  { m_model.load_state(state,overflow_state);  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) { m_model.store_state(state,overflow_state); }

  int adjoint_num_params(void) const { return material_adjoint_num_params(m_model); }

  const char* adjoint_param_name(int i) const { return material_adjoint_param_name(m_model,i); }

  void adjoint_add_stress_seed(Real* state, Real seed) { material_adjoint_add_stress_seed(m_model,state,seed); }

  void adjoint_objective_seed(Real* state) { material_adjoint_objective_seed(m_model,state); }

  Real adjoint_get_param_gradient(const Real* state, int i) const {
    return material_adjoint_get_param_gradient(m_model,state,i);
  }

  void adjoint_clear_step_seed(Real* state) { material_adjoint_clear_step_seed(m_model,state); }

  void adjoint_pullback_stress_to_strain(const Real* state, Real stress_seed, Real& strain_seed) const {
    material_adjoint_pullback_stress_to_strain(m_model,state,stress_seed,strain_seed);
  }

  void adjoint_add_direct_param_seed_from_stress(Real* state, Real stress_seed) {
    material_adjoint_add_direct_param_seed_from_stress(m_model,state,stress_seed);
  }

  Real get_state_variable(Real* state, std::string state_variable_name) { return m_model.get_state_variable(state,state_variable_name); }
  
  // Report the element death status of the current element
  bool is_dead(Real* state) { return m_model.is_dead(state); }
    
  // Return the number of nodes per element
  constexpr int num_nodes(void) { return 2; }

}; // Truss

#endif // TRUSS_H
