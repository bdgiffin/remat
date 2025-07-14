#ifndef UNIAXIAL_MATERIAL_H
#define UNIAXIAL_MATERIAL_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"

class UniaxialMaterial {
 protected:
  Real rho;   // Mass density
  Real area;  // Cross-sectional area
  Real E;     // Young's modulus
 public:

  // Empty constructor
  UniaxialMaterial(void) { }

  // Parameterized constructor
  UniaxialMaterial(Parameters& params) {
    // Get material density
    if (params.count("density") > 0) {
      rho = params["density"];
    } else {
      std::cout << "Missing uniaxial material parameter: density" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Get cross-sectional area
    if (params.count("area") > 0) {
      area = params["area"];
    } else {
      std::cout << "Missing uniaxial material parameter: area" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get primary elastic constants
    if (params.count("youngs_modulus") > 0) {
      E = params["youngs_modulus"];
    } else {
      std::cout << "Missing uniaxial material parameter: youngs_modulus" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return 1; }

  // Return the (mass per unit length) = (cross-sectional area) * (density)
  Real mass_per_unit_length(void) { return area*rho; }
    
  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0; // axial force
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

    // Compute and store the axial stress and force
    Real stress = E*strain;
    state[0] = stress*area;

    // Compute the strain energy per unit length
    psi = 0.5*state[0]*strain;
    
  } // update()

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(E/rho); }

}; /* UniaxialMaterial */

#endif /* UNIAXIAL_MATERIAL_H */
