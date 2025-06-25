#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "types.h"
#include <iostream>
#include <math.h>

template<class FixedV, class FixedU, class Ratio>
struct Integrator {
    
  // ================= INITIAL HALF-STEP VELOCITY UPDATE ================= //

  // Procedure to update first half-step velocities
  template<class Tv>
  void first_half_step_velocity_update(Real dt, Tv* v, Real* a, Real* alpha, int N) {
    
    // Loop over all degrees of freedom N
    for (int i=0; i<N; i++) {

      // 1) Compute the half-step velocity increment (cast as a FixedV number)
      FixedV dv = + 0.5*dt*a[i];

      // 2) Compute the damping factor for the current node (cast as a Rational number)
      Ratio phi = 1 - 0.5*dt*alpha[i];

      // 3) Update the dual velocity variable
      v[i] -= Dual<FixedV>(0,dv);

      // 4) Apply the damping factor to the primal/dual velocity pair
      v[i] *= phi;

      // 5) Update the primal velocity variable
      v[i] += Dual<FixedV>(dv,0);

    } // End loop over all DoFs
    
  } // initial_half_step_velocity_update()
    
  // ================== FINAL HALF-STEP VELOCITY UPDATE ================== //

  // Procedure to update second half-step velocities
  template<class Tv>
  void second_half_step_velocity_update(Real dt, Tv* v, Real* a, Real* alpha, int N) {

    // Loop over all degrees of freedom N
    for (int i=0; i<N; i++) {

      // 1) Compute the half-step velocity increment (cast as a FixedV number)
      FixedV dv = - 0.5*dt*a[i];

      // 2) Compute the damping factor for the current node (cast as a Rational number)
      Ratio phi = 1 + 0.5*dt*alpha[i];

      // 3) Update the dual velocity variable
      v[i] -= Dual<FixedV>(dv,0);

      // 4) Apply the damping factor to the primal/dual velocity pair
      v[i] /= phi;

      // 5) Update the primal velocity variable
      v[i] += Dual<FixedV>(0,dv);

    } // End loop over all DoFs
    
  } // final_half_step_velocity_update()
    
  // =================== WHOLE-STEP DISPLACEMENT UPDATE ================== //

  // Procedure to update whole-step displacements
  template<class Tv, class Tu>
  void whole_step_displacement_update(Real dt, Tv* v, Tu* u, int N) {
    
    // Loop over all degrees of freedom N
    for (int i=0; i<N; i++) {

      // 1) Compute the displacement increment (cast as a FixedU number) from the mid-step velocity
      FixedU du = v[i]*dt;
    
      // 2) Update the displacement with the displacement increment
      u[i] += Dual<FixedU>(+du,-du);

    } // End loop over all DoFs
    
  } // whole_step_displacement_update()
    
  // ===================================================================== //

private:
  
  // private constructor method for a time-integrator object
  //Integrator() { }
  
}; // Integrator

#endif // INTEGRATOR_H

