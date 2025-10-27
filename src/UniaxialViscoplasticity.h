#ifndef UNIAXIAL_VISCOPLASTICITY_H
#define UNIAXIAL_VISCOPLASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include <limits>
#include <vector>
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

  // Return the names of all fields
  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({ "axial_force",
                                      "axial_strain",
	                              "plastic_strain",
	                              "dual_plastic_strain",
	                              "overflow_counter",
	                              "equivalent_plastic_strain",
	                              "steps_since_element_death" });
  }

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
    Dual<FixedE> old_plastic_strain = plastic_strain;

    // Load the overflow counter from memory
    int overflow_counter = int(state[4]);

    // Load the equivalent plastic strain from memory
    FixedE eqps;
    load_from_Real(state[5],eqps);

    // Load the steps since element death
    int steps_since_death = int(state[6]);

    // Real epsilon_prev = state[1];
    if (dt > 0.0) {
      // epsilon_prev = strain;
      state[1] = strain;
    } 

    if (steps_since_death > 0) {
      if (dt >= 0.0) {
        steps_since_death++;
      } else {
    steps_since_death--;
      }
      std::cout << " since death = " << steps_since_death << std::endl;
    }

    Dual<FixedE> epsilon_prev_dual(FixedE(state[1]), FixedE(0.0));
    Dual<FixedE> elastic_strain = epsilon_prev_dual - old_plastic_strain;
    Real elastic_strain_primal = Real(elastic_strain.first);

    Real yield_strain = (E > 0.0) ? yield / E : 0.0;

    // Define a small tolerance in the fixed-point representation
    Dual<FixedE> n_times_epsilon_tol_fixed(FixedE(0.0), FixedE(0.0));
    
    // Check for yielding
    if ((std::abs(elastic_strain_primal) >= yield_strain) && (steps_since_death == 0)) {
      Real n_sign = 0.0;
      if (elastic_strain_primal > 0.0) {
        n_times_epsilon_tol_fixed.first = FixedE(1e-6);
        n_sign = 1.0;
      } else if (elastic_strain_primal < 0.0) {
        n_sign = -1.0;
        n_times_epsilon_tol_fixed.first = FixedE(-1e-6);
      }

      Real A = std::exp(-E*std::abs(dt)/eta);
      Ratio A_rat(A);

      // Real offset_real = n_sign * (yield_strain + epsilon_tol_real);
      Real offset_real = n_sign * (yield_strain) ;
      Dual<FixedE> offset(FixedE(offset_real), FixedE(0.0));
      offset = offset + n_times_epsilon_tol_fixed;
      Dual<FixedE> delta_elastic = elastic_strain - offset;

      if (dt > 0.0) {
        delta_elastic = delta_elastic * A_rat;
      } else {
        delta_elastic = delta_elastic / A_rat;
      }

      elastic_strain = delta_elastic + offset;
      plastic_strain = epsilon_prev_dual - elastic_strain;

      std::cout << "yielding" << std::endl;
      if (dt >= 0.0) {
        FixedE delta_ps = plastic_strain.first - old_plastic_strain.first;
	      if (delta_ps < FixedE(0.0)) {
          eqps = eqps - delta_ps;
        } else {
          eqps = eqps + delta_ps;
        }
        overflow_counter++;
        if (Real(eqps) > epsf) { steps_since_death = 1; }
      } else {
        overflow_counter--;
        FixedE delta_ps = old_plastic_strain.first - plastic_strain.first;
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
    plastic_strain = epsilon_prev_dual - elastic_strain;

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
  void load_state(Real* state, std::vector<Real>& overflow_state)  {
    if ((int(state[4]) == 0) && (overflow_state.size() > 0)) {
      state[4] = Real(mat_overflow_limit);
      state[3] = overflow_state.back(); overflow_state.pop_back();
    }
  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    if (int(state[4]) == mat_overflow_limit) {
      state[4] = Real(0);
      overflow_state.push_back(state[3]);
      state[3] = smallest_value(state[3]);
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

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    FixedE temp;
    field_data[0] = state[0]; // axial_force
    field_data[1] = state[1]; // axial_strain
    load_from_Real(state[2],temp); field_data[2] = Real(temp); // plastic_strain
    load_from_Real(state[3],temp); field_data[3] = Real(temp); // dual_plastic_strain
    field_data[4] = state[4]; // overflow_counter
    load_from_Real(state[5],temp); field_data[5] = Real(temp); // equivalent_plastic_strain
    field_data[6] = state[6]; // steps_since_element_death
  }

  // Report the element death status of the current material point
  bool is_dead(Real* state) { return (int(state[6]) > 0); }

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(E/rho); }

}; /* UniaxialViscoplasticity */

#endif /* UNIAXIAL_VISCOPLASTICITY_H */
