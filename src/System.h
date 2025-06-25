#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.h"
#include <vector>
#include <iostream>
#include <math.h>

template<class FixedV, class FixedU>
struct System {
  int N;                        // The total number of system degrees of freedom
  std::vector<Dual<FixedU> > u; // The (primal/dual) system displacement degrees of freedom
  std::vector<Dual<FixedV> > v; // The (primal/dual) system velocity degrees of freedom
  std::vector<Real>          a; // The (primal) system accelerations for each DoF
  std::vector<Real>      alpha; // The system damping factors for each DoF
  Integrator<FixedV,FixedU,Rational> leapfrog; // (Bit-reversible) leapfrog time integrator
    
  // ===================================================================== //

  // Procedure to initialize the system state
  void update_accelerations(Real time, Real dt) {

    // Initialize forces and masses

    // Sum external forces

    // Loop over all elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal displacements for each element (cast as Reals)

      // Compute internal force in the current element, and sum global forces and masses

    } // End loop over all elements

    // Compute updated accelerations and damping factors for each DoF:

    
    // Loop over all nodes
    for (int i=0; i<Nnodes; i++) {

      // 1) Cast the primal displacement variables from FixedU numbers to floating point
      // (copy from global array into local array for each element)
      Real ur = u;
    
      // 2) Evaluate the residual forces at the current time state, and assemble nodal masses
      update_accelerations(time,dt);

      // 3) Divide nodal forces by nodal masses to obtain the predictor nodal accelerations
      a[i] = f[i]/m[i];

    } // End loop over all nodes
    
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
    
    // Compute initial masses and residual forces at the current (initial) time
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

