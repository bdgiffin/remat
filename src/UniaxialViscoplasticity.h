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
  Real tau;   // Relaxation time
  Real yield; // Yield stress
  Real epsf;  // Equivalent plastic strian at failure
  int  mat_overflow_limit = std::numeric_limits<int>::max();
  FixedE yield_strain; // Yield strain
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

    // Verify that the Young's modulus is a strictly positive value:
    if (E <= 0.0) {
      std::cout << "Invalid uniaxial viscoplasticity parameter specified: (youngs_modulus / truss_youngs_modulus) must be a strictly positive value" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get viscosity
    if (params.count("viscosity") > 0) {
      eta = params["viscosity"];
    } else {
      eta = std::numeric_limits<Real>::max();
    }

    // Compute the relaxation time
    tau = eta/E;

    // Get yield_stress
    if (params.count("yield_stress") > 0) {
      yield = params["yield_stress"];
    } else {
      yield = std::numeric_limits<Real>::max();
    }

    // Compute the initial yield strain
    yield_strain = FixedE(yield/E);

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
      std::cout << "ERROR: non-positive stretch ratio" << std::endl;
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
    
    // Increment/decrement the "steps since element death" counter
    if (steps_since_death > 0) {
      if (dt >= 0.0) {
        steps_since_death++;
      } else {
	steps_since_death--;
      }
    }

    // Conditionally update the strain history variable
    if (dt > 0.0) { state[1] = strain; } // strain_prev <- strain
    Dual<FixedE> total_strain(FixedE(state[1]), FixedE(0.0));

    // Compute the trial elastic strain = (total strain) - (plastic strain)
    Dual<FixedE> elastic_strain = total_strain - plastic_strain;
    
    // Check for yielding
    if ((abs(elastic_strain.first) >= yield_strain) && (steps_since_death == 0)) {

      // Determine the flow direction
      bool n_positive = (elastic_strain.first >= FixedE(0.0));

      // Compute the amplification factor
      Ratio A_rat(std::exp(-std::abs(dt)/tau));

      // Offset the trial elastic strain by the (signed) yield strain
      if (n_positive) {
	elastic_strain.first = elastic_strain.first - yield_strain;
      } else {
	elastic_strain.first = elastic_strain.first + yield_strain;
      }

      // Multiply the shifted elastic strain by the amplification factor
      FixedE old_elastic_strain = elastic_strain.first;
      if (dt > 0.0) { // forward-in-time
        elastic_strain = elastic_strain * A_rat;
      } else {        // backward-in-time
        elastic_strain = elastic_strain / A_rat;
      }

      // Determine the change in plastic strain (negative of the change in elastic strain)
      // Delta Ep = (old Ee) - (new Ee) = (new Ep) - (old Ep)
      // FixedE delta_ps = old_elastic_strain - elastic_strain.first;

      // Update the equivalent plastic strain
      FixedE delta_eqps = abs(old_elastic_strain - elastic_strain.first);
      if (dt >= 0.0) {
	eqps = eqps + delta_eqps;
        overflow_counter++;
        if (Real(eqps) > epsf) { steps_since_death = 1; }
      } else {
        overflow_counter--;
	eqps = eqps - delta_eqps;
      }

      // Reverse-offset the trial elastic strain by the (signed) yield strain
      if (n_positive) {
	elastic_strain.first = elastic_strain.first + yield_strain;
      } else {
	elastic_strain.first = elastic_strain.first - yield_strain;
      }
      
    } // check for yielding

    // Compute the updated plastic strain = (total strain) - (elastic strain)
    plastic_strain = total_strain - elastic_strain;

    // Conditionally update the strain history variable
    if (dt < 0.0) { state[1] = strain; } // strain_prev <- strain

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
      state[3] = Real(0.0);
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
