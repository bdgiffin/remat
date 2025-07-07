#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.h"
#include <vector>
#include <algorithm> // For std::fill
#include <iostream>
#include <math.h>

template<class FixedV, class FixedU>
struct System {
  int Nnodes;                   // The total number of nodes
  int Ndofs_per_node;           // The total number of degrees of freedom per node
  int Ndofs;                    // The total number of system degrees of freedom = Nnodes * Ndofs_per_node
  int Nelems;                   // The total number of elements
  int Nnodes_per_elem;          // The total number of nodes per element
  std::vector<int>     connect; // The nodal connectivity array for all elements
  std::vector<Real>          x; // The (primal) system initial position degrees of freedom
  std::vector<Dual<FixedU> > u; // The (primal/dual) system displacement degrees of freedom
  std::vector<Dual<FixedV> > v; // The (primal/dual) system velocity degrees of freedom
  std::vector<Real>          m; // The (primal) system masses for each DoF
  std::vector<Real>          f; // The (primal) system forces for each DoF
  std::vector<Real>          a; // The (primal) system accelerations for each DoF
  std::vector<Real>      alpha; // The system damping factors for each DoF
  Integrator<FixedV,FixedU,Rational> leapfrog; // (Bit-reversible) leapfrog time integrator
    
  // ===================================================================== //

  // Procedure to initialize the system
  void initialize(Real* new_coordinates, int new_Nnodes, int new_Ndofs_per_node,
		  int* new_connectivity, int new_Nelems, int new_Nnodes_per_elem) {

    // Initialize all mesh totals
    Nnodes          = new_Nnodes;
    Ndofs_per_node  = new_Ndofs_per_node;
    Ndofs           = Nnodes*Ndofs_per_node;
    Nelems          = new_Nelems;
    Nnodes_per_elem = new_Nnodes_per_elem;
    const int Nelem_dofs = Nelems*Nnodes_per_elem;

    // Initialize the dimensions of all arrays
    connect.resize(Nelem_dofs);
    x.resize(Ndofs);
    u.resize(Ndofs);
    v.resize(Ndofs);
    m.resize(Ndofs);
    f.resize(Ndofs);
    a.resize(Ndofs);
    alpha.resize(Ndofs);

    // Initialize data for all elements
    for (int i=0; i<Nelem_dofs; i++) {
      connect[i] = new_connectivity[i];
    }

    // Initialize data for all DoFs
    for (int i=0; i<Ndofs; i++) {
      x[i] = new_coordinates[i];
      u[i] = Dual<FixedU>(0,0);
      v[i] = Dual<FixedV>(0,0);
      m[i] = 0;
      f[i] = 0;
      a[i] = 0;
      alpha[i] = 0;
    }
    
  } // initialize()
  
  // ===================================================================== //

  // Procedure to update the accelerations
  void update_accelerations(Real time, Real dt) {

    // Initialize forces and masses
    std::fill(m.begin(), m.end(), 0);
    std::fill(f.begin(), f.end(), 0); 

    // Sum external forces
    // FIXME: add external forces

    // Loop over all elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal displacements for each element (cast as Reals)
      const int Ndofs_per_elem = Nnodes_per_elem*Ndofs_per_node;
      Real xe[Ndofs_per_elem];
      Real ue[Ndofs_per_elem];
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
	  ue[Ndofs_per_node*j+i] = u[Ndofs_per_node*jnode_id+i];
	}
      }

      // Compute internal force in the current element
      Real me[Ndofs_per_elem];
      Real fe[Ndofs_per_elem];
      compute_element_mass_and_forces(xe,ue,me,fe);

      // Scatter mass and forces to the nodes
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  m[Ndofs_per_node*jnode_id+i] += me[Ndofs_per_node*j+i];
	  f[Ndofs_per_node*jnode_id+i] += fe[Ndofs_per_node*j+i];
	}
      }

    } // End loop over all elements

    // Compute updated accelerations and damping factors for each DoF:
    for (int i=0; i<Ndofs; i++) {

      // Divide nodal forces by nodal masses to obtain the (predictor) nodal accelerations
      a[i] = f[i]/m[i];

      // Determine the (mass-proportional) damping factor "alpha" in terms of a target frequency and damping ratio
      Real frequency = 1.0;
      Real damping_ratio = 0.1;
      alpha[i] = 2.0*frequency*damping_ratio;

    } // End loop over all DoFs
    
  } // update_accelerations()
  
  // ===================================================================== //

  // Procedure to initialize the system state
  void initialize_state(Real time) {
    
    // Update accelerations and damping factors for each DoF
    update_accelerations(time,0);
    
  } // initialize_state()
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  void update_state(Real time, Real dt) {

    // Update velocities to the half-step
    leapfrog.first_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),N);

    // Update displacement to the next whole-step
    leapfrog.whole_step_displacement_update(dt,v.data(),u.data(),N);
    
    // Update masses, residual forces, and accelerations at the whole-step
    update_accelerations(time,dt);
    
    // Update velocities to the whole-step
    leapfrog.second_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),N);
    
  } // update_state()
  
  // ===================================================================== //

private:
  
  // private constructor method for a system object
  System() { }
  
}; // System

#endif // SYSTEM_H

