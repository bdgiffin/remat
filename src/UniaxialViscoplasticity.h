#ifndef UNIAXIAL_VISCOPLASTICITY_H
#define UNIAXIAL_VISCOPLASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include <limits>
#include "types.h"
#include "Dual.h"
#include "Parameters.h"

template<class FixedE, class Ratio>
class UniaxialViscoplasticity {
 protected:
  Real rho;   // Mass density
  Real area;  // Cross-sectional area
  Real E;     // Young's modulus
  Real eta;   // Viscosity
  Real yield; // Yield stress
  Real epsf;  // Equivalent plastic strian at failure
  int  mat_overflow_limit = std::numeric_limits<int>::max();
 public:

  // Empty constructor
  UniaxialViscoplasticity(void) { }

  // Parameterized constructor
  UniaxialViscoplasticity(Parameters& params) {
    // Get material density
    if  (params.count("truss_density") > 0) {
      rho = params["truss_density"];
    } else if (params.count("density") > 0) {
      rho = params["density"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: density / truss_density" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Get cross-sectional area
    if (params.count("area") > 0) {
      area = params["area"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: area" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get primary elastic constants
    if  (params.count("truss_youngs_modulus") > 0) {
      E = params["truss_youngs_modulus"];
    } else if (params.count("youngs_modulus") > 0) {
      E = params["youngs_modulus"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: youngs_modulus / truss_youngs_modulus" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get viscosity
    if (params.count("viscosity") > 0) {
      eta = params["viscosity"];
    } else {
      eta = std::numeric_limits<Real>::max();
    }

    // Get yield_stress
    if (params.count("yield_stress") > 0) {
      yield = params["yield_stress"];
    } else {
      yield = std::numeric_limits<Real>::max();
    }

    // Get eps_fail (equivalent plastic strain at failure)
    if (params.count("eps_fail") > 0) {
      epsf = params["eps_fail"];
    } else {
      epsf = std::numeric_limits<Real>::max();
    }

    // Get overflow limit for material history variables
    if (params.count("mat_overflow_limit") > 0) { mat_overflow_limit = int(params["mat_overflow_limit"]); }
  }

  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return 7; }

  // Return the (mass per unit length) = (cross-sectional area) * (density)
  Real mass_per_unit_length(void) { return area*rho; }
  
  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0; // axial force
    state[1] = 0.0; // axial strain
    save_as_Real(FixedE(0.0),state[2]); // plastic strain
    save_as_Real(FixedE(0.0),state[3]); // dual plastic strain
    state[4] = Real(0); // overflow counter
    save_as_Real(FixedE(0.0),state[5]); // equivalent plastic strain
    state[6] = Real(0); // steps since element death
  } // initialize()
  
  // Update the material state using the current stretch ratio
  void update(Real lambda, Real &psi, Real* state, Real dt) {
    // Simple linear elastic model model:

    // Check for element inversion
    if (lambda <= 0.0) {
      std::cout << "ERROR: negative stretch ratio" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute the axial (small) strain
    Real strain = lambda - 1.0;

    // Load the plastic strain and its dual from memory
    FixedE p1,p2;
    load_from_Real(state[2],p1);
    load_from_Real(state[3],p2);
    Dual<FixedE> plastic_strain(p1,p2);

    // Load the overflow counter from memory
    int overflow_counter = int(state[4]);

    // Load the equivalent plastic strain from memory
    FixedE eqps;
    load_from_Real(state[5],eqps);

    // Load the steps since element death
    int steps_since_death = int(state[6]);

    // Compute the scaling factor
    Real phi = std::exp(-E*std::abs(dt)/eta);
    Ratio rat(phi);

    // Conditionally update the strain history variable
    if (dt >= 0.0) state[1] = strain;

    // Determine the flow direction and the reference plastic strain
    Real n = std::copysign(1.0,state[1] - Real(plastic_strain.first));
    Real ep0 = state[1] - n*yield/E;

    // Update the number of steps since element death occurred
    if (steps_since_death > 0) {
      if (dt >= 0.0) {
        steps_since_death++;
      } else {
	steps_since_death--;
      }
      std::cout << " since death = " << steps_since_death << std::endl;
    }

    // Check for yielding
    if ((std::abs(E*(state[1] - Real(plastic_strain.first))) > yield) && (steps_since_death == 0)) {
      // Update the plastic strain
      Dual<FixedE> dep(ep0*(1.0 - phi), 0.0);
      std::cout << "yielding" << std::endl;
      if (dt >= 0.0) {
	FixedE old_ps = plastic_strain.first;
	plastic_strain = plastic_strain*rat;
	plastic_strain = plastic_strain + dep;
	FixedE delta_ps = plastic_strain.first - old_ps;
	if (delta_ps < FixedE(0.0)) {
	  eqps = eqps - delta_ps;
	} else {
	  eqps = eqps + delta_ps;
	}
	overflow_counter++;
	if (Real(eqps) > epsf) { steps_since_death = 1; } // the element has died
      } else {
	overflow_counter--;
	FixedE new_ps = plastic_strain.first;
	plastic_strain = plastic_strain - dep;
	plastic_strain = plastic_strain/rat;
	FixedE delta_ps = new_ps - plastic_strain.first;
	if (delta_ps < FixedE(0.0)) {
	  eqps = eqps + delta_ps;
	} else {
	  eqps = eqps - delta_ps;
	}
      }
      //std::cout << " growth fac = " << 1.0/phi << std::endl;
      //std::cout << "  primal ps = " << Real(plastic_strain.first)  << std::endl;
      //std::cout << "    dual ps = " << Real(plastic_strain.second) << std::endl;
      //std::cout << "    counter = " << overflow_counter << std::endl;
      std::cout << "    eqps = " << Real(eqps) << std::endl;
    }

    // Conditionally update the strain history variable
    if (dt < 0.0) state[1] = strain;
    
    // Compute and store the axial stress and force
    Real stress = (steps_since_death > 0) ? 0.0 : E*(state[1] - Real(plastic_strain.first));
    state[0] = stress*area;

    // Save the plastic strain and its dual to memory
    save_as_Real(plastic_strain.first,  state[2]);
    save_as_Real(plastic_strain.second, state[3]);

    // Save the overflow counter to memory
    state[4] = Real(overflow_counter);

    // Save the equivalent plastic strain to memory
    save_as_Real(eqps, state[5]);

    // Save the steps since element death
    state[6] = Real(steps_since_death);

    // Compute the strain energy per unit length
    psi = (steps_since_death > 0) ? 0.0 : 0.5*state[0]*(state[1] - Real(plastic_strain.first));
    
  } // update()

  // Conditionally load material history parameters from memory
  bool load_state(Real* state, Real* overflow_state)  {
    if (int(state[4]) == 0) {
      state[4] = Real(mat_overflow_limit);
      state[3] = *overflow_state;
      return true;
    } else {
      return false;
    }
  }

  // Conditionally store material history parameters in memory
  bool store_state(Real* state, Real* overflow_state) {
    if (int(state[4]) == mat_overflow_limit) {
      state[4] = Real(0);
      *overflow_state = state[3];
      state[3] = smallest_value(state[3]);
      return true;
    } else {
      return false;
    }
  }
  
  Real get_state_variable(Real* state, std::string state_variable_name) {
    if (state_variable_name == "eqps") {
      FixedE eqps;
      load_from_Real(state[5],eqps);
      return Real(eqps);
    } else {
      return Real(0.0);
    }
  }

  // Report the element death status of the current material point
  bool is_dead(Real* state) { return (int(state[6]) > 0); }

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(E/rho); }

}; /* UniaxialViscoplasticity */

#endif /* UNIAXIAL_VISCOPLASTICITY_H */
