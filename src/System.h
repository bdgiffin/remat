#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.h"
#include "Dual.h"
#include "Integrator.h"
#include "Parameters.h"
#include <vector>
#include <algorithm> // For std::fill
#include <iostream>
#include <math.h>

template<class Element_T, class FixedV, class FixedU, class Ratio>
struct System {

  // Data members:

  // Time step ID
  int m_time_step;

  // Analysis time
  Real m_time;

  // Stable time step size
  Real m_dt_scale_factor;
  Real m_stable_dt;

  // Global body forces
  Real m_bx = 0.0;
  Real m_by = 0.0;

  // Mass-proportional damping factor
  Real m_alpha = 0.0;

  // Contact stiffness parameter
  Real m_contact_stiffness = 0.0;

  // Initial velocity
  Real m_vx0 = 0.0;
  Real m_vy0 = 0.0;

  // Mesh object and data
  //Mesh m_mesh;

  // Boundary conditions
  //std::vector<BoundaryCondition> m_bcs;
  //std::vector<NodalForce> m_forces;
  
  int Nnodes;                   // The total number of nodes
  int Ndofs_per_node;           // The total number of degrees of freedom per node
  int Ndofs;                    // The total number of system degrees of freedom = Nnodes * Ndofs_per_node
  int Nelems;                   // The total number of elements
  int Nnodes_per_elem;          // The total number of nodes per element
  std::vector<int>     connect; // The nodal connectivity array for all elements
  std::vector<bool>     fixity; // The fixity (fixed == true) assigned to each nodal degree of freedom
  std::vector<Real>          x; // The (primal) system initial position degrees of freedom
  std::vector<Real>         xt; // The (primal) system current position degrees of freedom
  std::vector<Dual<FixedU> > u; // The (primal/dual) system displacement degrees of freedom
  std::vector<Dual<FixedV> > v; // The (primal/dual) system velocity degrees of freedom
  std::vector<Real>          m; // The (primal) system masses for each DoF
  std::vector<Real>          f; // The (primal) system forces for each DoF
  std::vector<Dual<Real> >   a; // The (primal) system accelerations for each DoF
  std::vector<Real>      alpha; // The system damping factors for each DoF
  std::vector<Real>      state; // Element state variable data
  
  Integrator<Ratio> m_integrator;                         // (Bit-reversible) leapfrog time integrator
  Element_T m_element;                                    // Element class
  std::vector<ContactInteraction> m_contact_interactions; // List of penalty-based contact interactions
    
  // ===================================================================== //
  
  // default constructor method for a system object
  System() { }
    
  // ===================================================================== //

  // Procedure to initialize the system
  void initialize(Real* new_coordinates, Real* new_velocities, bool* new_fixity, int new_Nnodes, int new_Ndofs_per_node,
		  int* new_connectivity, int new_Nelems, int new_Nnodes_per_elem,
		  Parameters& params) {
    
    std::cout << "| ====================== INITIALIZING ====================== |" << std::endl;

    std::cout << "Initializing state data..." << std::endl;

    // Initialize global parameters
    m_time_step = 0;
    m_time = 0.0;
    if (params.count("body_force_x") > 0) m_bx = params["body_force_x"];
    if (params.count("body_force_y") > 0) m_by = params["body_force_y"];
    if (params.count("initial_velocity_x") > 0) m_vx0 = params["initial_velocity_x"];
    if (params.count("initial_velocity_y") > 0) m_vy0 = params["initial_velocity_y"];
    if (params.count("dt_scale_factor") > 0) m_dt_scale_factor = params["dt_scale_factor"];
    if (params.count("mass_damping_factor") > 0) m_alpha = params["mass_damping_factor"];
    if (params.count("contact_stiffness") > 0) m_contact_stiffness = params["contact_stiffness"];

    // Initialize the element/material object
    m_element = Element_T(params);

    // Initialize all mesh totals
    Nnodes          = new_Nnodes;
    Ndofs_per_node  = new_Ndofs_per_node;
    Ndofs           = Nnodes*Ndofs_per_node;
    Nelems          = new_Nelems;
    Nnodes_per_elem = new_Nnodes_per_elem;
    const int Nelem_dofs = Nelems*Nnodes_per_elem;
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    const int Nstate = Nelems*Nstate_vars_per_elem;

    // Initialize the dimensions of all arrays
    connect.resize(Nelem_dofs);
    fixity.resize(Ndofs);
    x.resize(Ndofs);
    xt.resize(Ndofs);
    u.resize(Ndofs,Dual<FixedU>(0.0,0.0));
    v.resize(Ndofs,Dual<FixedV>(0.0,0.0));
    m.resize(Ndofs);
    f.resize(Ndofs);
    a.resize(Ndofs,Dual<Real>(0.0,0.0));
    alpha.resize(Ndofs);
    state.resize(Nstate);

    // Initialize data for all DoFs
    for (int i=0; i<Ndofs; i++) {
      fixity[i] = new_fixity[i];
      x[i] = new_coordinates[i];
      xt[i] = x[i];
      u[i] = Dual<FixedU>(0.0,0.0);
      v[i] = Dual<FixedV>(new_velocities[i],0.0);
      m[i] = 0.0;
      f[i] = 0.0;
      a[i] = Dual<Real>(0.0,0.0);
      alpha[i] = 0.0;
    }
    
    // Assign constant initial velocity (if defined)
    for (int i = 0; i < Nnodes; i++) {
      if (m_vx0 != 0.0) v[2*i+0] = Dual<FixedV>(m_vx0,0.0);
      if (m_vy0 != 0.0) v[2*i+1] = Dual<FixedV>(m_vy0,0.0);
    }

    // Initialize connectivity data for all elements
    for (int i=0; i<Nelem_dofs; i++) {
      connect[i] = new_connectivity[i];
    }

    // Initialize state variable data for all elements
    for (int e=0; e<Nelems; e++) {
      m_element.initialize(&state[Nstate_vars_per_elem*e]);
    }
    
    std::cout << "| ========================================================== |" << std::endl;
    
  } // initialize()
  
  // ===================================================================== //

  // Procedure to initialize a new contact interaction
  void initialize_contact(int* node_ids, int* segment_connectivity, int new_Nnodes, int new_Nsegments, Parameters& params) {
    
    // Define a new contact interaction, and initialize it
    m_contact_interactions.push_back(ContactInteraction());
    m_contact_interactions.back().initialize(node_ids,segment_connectivity,new_Nnodes,new_Nsegments,params);
    
  } // initialize_contact()
  
  // ===================================================================== //

  // Procedure to initialize the system state at the indicated time
  void initialize_state(void) {

    // Set the initial time to zero
    m_time = 0.0;
    
    // Update accelerations and damping factors for each DoF
    update_accelerations(0.0);
    
  } // initialize_state()
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  Real update_state(Real dt) {

    // Update velocities to the half-step
    m_integrator.first_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),Ndofs);

    // Update displacement to the next whole-step
    m_integrator.whole_step_displacement_update(dt,v.data(),u.data(),Ndofs);
    
    // Update masses, residual forces, and accelerations at the whole-step
    update_accelerations(dt);
    
    // Update velocities to the whole-step
    m_integrator.second_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),Ndofs);

    // Update time step ID
    if (dt > 0.0) {
      m_time_step++;
    } else if (dt < 0.0) {
      m_time_step--;
    }
    std::cout << "Time step: " << m_time_step << " at time: " << m_time << std::endl;

    // Return the updated analysis time
    return m_time;
    
  } // update_state()
  
  // ===================================================================== //

  // Procedure to get the current system state data
  double get_field_data(double *ux, double *uy,
	 	        double *vx, double *vy,
		        double *fx, double *fy,
			double *dual_ux, double *dual_uy,
	 	        double *dual_vx, double *dual_vy,
		        double *sxx, double *syy, double *sxy, double *pressure) {

    // Copy nodal state data
    for (int i=0; i<Nnodes; i++) {
      ux[i]      = u[2*i+0].first;
      uy[i]      = u[2*i+1].first;
      vx[i]      = v[2*i+0].first;
      vy[i]      = v[2*i+1].first;
      fx[i]      = f[2*i+0];
      fy[i]      = f[2*i+1];
      dual_ux[i] = u[2*i+0].second;
      dual_uy[i] = u[2*i+1].second;
      dual_vx[i] = v[2*i+0].second;
      dual_vy[i] = v[2*i+1].second;
    }

    // Copy element state data
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    for (int e=0; e<Nelems; e++) {
      sxx[e]      = state[Nstate_vars_per_elem*e+0];
      syy[e]      = state[Nstate_vars_per_elem*e+1];
      sxy[e]      = state[Nstate_vars_per_elem*e+2];
      pressure[e] = state[Nstate_vars_per_elem*e+3];
    }

    // Return the current analysis time
    return m_time;
    
  } // update_state()
  
  // ===================================================================== //
private:
  // ===================================================================== //

  // Procedure to update the accelerations
  void update_accelerations(Real dt) {

    // Update the current analysis time
    m_time = m_time + dt;

    // Zero-initialize forces and masses
    std::fill(m.begin(), m.end(), 0.0);
    std::fill(f.begin(), f.end(), 0.0);

    // Update the current deformed nodal coordinates
    for (int i = 0; i < Nnodes; i++) {
      xt[2*i+0] = x[2*i+0] + u[2*i+0].first;
      xt[2*i+1] = x[2*i+1] + u[2*i+1].first;
    }
    
    const int Nstate_vars_per_elem = m_element.num_state_vars();

    // Loop over all elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal displacements for each element (cast as Reals)
      const int Ndofs_per_elem = Nnodes_per_elem*Ndofs_per_node;
      Real xe[8];
      Real ue[8];
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
	  ue[Ndofs_per_node*j+i] = u[Ndofs_per_node*jnode_id+i].first;
	}
      }

      // Compute internal force in the current element
      Real me[8] = { 0.0 };
      Real fe[8] = { 0.0 };
      m_element.update(xe,ue,me,fe,&state[Nstate_vars_per_elem*e],dt);

      // Scatter mass and forces to the nodes
      // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  m[Ndofs_per_node*jnode_id+i] += me[Ndofs_per_node*j+i];
	  f[Ndofs_per_node*jnode_id+i] -= fe[Ndofs_per_node*j+i];
	}
      }

    } // End loop over all elements

    // Sum nodal forces from external nodal point loads (gravity, etc.)
    for (int i = 0; i < Nnodes; i++) {
      // Global body forces
      f[2*i+0] += m_bx*m[2*i+0];
      f[2*i+1] += m_by*m[2*i+1];
    }

    // Sum nodal forces due to contact interactions
    for (auto& contact : m_contact_interactions) contact.update_contact_forces(xt,f,dt);

    // Contact interaction with rigid wall
    for (int i = 0; i < Nnodes; i++) {
      if (x[2*i+1] + u[2*i+1].first < 0.0) {
	Real contact_stiffness = 1.0e1;
	Real du = - (x[2*i+1] + u[2*i+1].first);
	f[2*i+1] += m_contact_stiffness*du;
      }
    }

    // Compute updated accelerations and damping factors for each DoF:
    for (int i=0; i<Ndofs; i++) {

      // Divide nodal forces by nodal masses to obtain the (predictor) nodal accelerations
      a[i].first = f[i]/m[i];

      // Determine fictitious dual nodal accelerations (directly proportional to the dual displacements)
      a[i].second = 0.0*u[i].second;

      // Determine the (mass-proportional) damping factor "alpha" in terms of a target frequency and damping ratio
      //Real frequency = 1.0;
      //Real damping_ratio = 0.1;
      //alpha[i] = 2.0*frequency*damping_ratio;
      alpha[i] = m_alpha;

      // Impose fixed boundary conditions
      if (fixity[i]) v[i].first = 0.0;

    } // End loop over all DoFs

    // Report the largest dual variable value (indicates an integer overflow ...)
    //Integer max_val;
    //for (int i = 0; i < Ndofs; i++) {
    //  max_val = std::max(max_val,std::abs(v[i].second.mantissa));
    //}
    //if (Real(max_val) > 0.99*std::numeric_limits<Integer>::max()) std::cout << "max val = " << max_val << " > " << 0.99*std::numeric_limits<Integer>::max() << std::endl;
    
  } // update_accelerations()
  
  // ===================================================================== //
  
}; // System

#endif // SYSTEM_H

