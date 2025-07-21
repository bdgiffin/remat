#ifndef UNIAXIAL_VISCOPLASTICITY_H
#define UNIAXIAL_VISCOPLASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
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
 public:

  // Empty constructor
  UniaxialViscoplasticity(void) { }

  // Parameterized constructor
  UniaxialViscoplasticity(Parameters& params) {
    // Get material density
    if (params.count("density") > 0) {
      rho = params["density"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: density" << std::endl;
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
    if (params.count("youngs_modulus") > 0) {
      E = params["youngs_modulus"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: youngs_modulus" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get viscosity
    if (params.count("viscosity") > 0) {
      eta = params["viscosity"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: viscosity" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get yield_stress
    if (params.count("yield_stress") > 0) {
      yield = params["yield_stress"];
    } else {
      std::cout << "Missing uniaxial viscoplasticity parameter: yield_stress" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return 4; }

  // Return the (mass per unit length) = (cross-sectional area) * (density)
  Real mass_per_unit_length(void) { return area*rho; }
    
  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0; // axial force
    state[1] = 0.0; // axial strain
    save_as_Real(FixedE(0.0),state[2]); // plastic strain
    save_as_Real(FixedE(0.0),state[3]); // dual plastic strain
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

    // Compute the scaling factor
    Real phi = std::exp(-E*std::abs(dt)/eta);
    Ratio rat(phi);

    // Conditionally update the strain history variable
    if (dt >= 0.0) state[1] = strain;

    // Determine the flow direction and the reference plastic strain
    Real n = std::copysign(1.0,state[1] - plastic_strain.first);
    Real ep0 = state[1] - n*yield/E;

    // Check for yielding
    if (std::abs(E*(state[1] - plastic_strain.first)) > yield) {
      // Update the plastic strain
      Dual<FixedE> dep(ep0*(1.0 - phi), 0.0);
      std::cout << "yielding" << std::endl;
      if (dt >= 0.0) {
	plastic_strain = plastic_strain*rat;
	plastic_strain = plastic_strain + dep;
      } else {
	plastic_strain = plastic_strain - dep;
	plastic_strain = plastic_strain/rat;
      }
      std::cout << " growth fac = " << 1.0/phi << std::endl;
      std::cout << "  primal ps = " << plastic_strain.first  << std::endl;
      std::cout << "    dual ps = " << plastic_strain.second << std::endl;
    }

    // Conditionally update the strain history variable
    if (dt < 0.0) state[1] = strain;
    
    // Compute and store the axial stress and force
    Real stress = E*(state[1] - plastic_strain.first);
    state[0] = stress*area;

    // Save the plastic strain and its dual to memory
    save_as_Real(plastic_strain.first,  state[2]);
    save_as_Real(plastic_strain.second, state[3]);

    // Compute the strain energy per unit length
    psi = 0.5*state[0]*(state[1] - plastic_strain.first);
    
  } // update()

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(E/rho); }

}; /* UniaxialViscoplasticity */

#endif /* UNIAXIAL_VISCOPLASTICITY_H */
