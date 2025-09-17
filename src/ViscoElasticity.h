#ifndef VISCOELASTICITY_H
#define VISCOELASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"

class ViscoElasticity {
 protected:
  Real rho;   // Mass density
  Real E;     // Young's modulus
  Real nu;    // Poisson's ratio
  Real mu;    // Shear modulus
  Real mu2;   // Twice the shear modulus
  Real lam;   // Lame parameter
  Real kappa; // Bulk modulus
  Real pmod;  // P-wave modulus
// Viscous parameters
  Real tau;    // relaxation time
  Real eta;    // viscosity (if provided then tau = eta/mu would be computed)
 public:

  // Empty constructor
  ViscoElasticity(void) { }

  // Parameterized constructor
  ViscoElasticity(Parameters& params) {
    // Get material density
    if (params.count("density") > 0) {
      rho = params["density"];
    } else {
      std::cout << "Missing material parameter: density" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get primary elastic constants
    if (params.count("youngs_modulus") > 0) {
      E = params["youngs_modulus"];
    } else {
      std::cout << "Missing material parameter: youngs_modulus" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (params.count("poissons_ratio") > 0) {
      nu = params["poissons_ratio"];
    } else {
      std::cout << "Missing material parameter: poissons_ratio" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute derived elastic constants
    mu2   = E/(1.0+nu);
    mu    = 0.5*mu2;
    lam   = mu2*nu/(1.0-2.0*nu);
    kappa = lam+0.5*mu2;
    pmod  = lam+mu2;

    // Visco parameters
    // Get relatxation_time if provided
    // Otherwise, relatxation_time = viscosity / G
    if (params.count("relaxation_time") > 0) {
      tau = params["relaxation_time"];
    } else {
      if (params.count("viscosity") > 0) {
        eta = params["viscosity"];
        tau = eta / mu;
      } else {      
        std::cout << "Missing visco parameters: relaxation_time" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (tau<=0.0) {
      std::cout << "ViscoElasticity: need a positive relaxation_time" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
  // Return the number of state variables for allocation purposes
 int num_state_vars(void) { return 12; }

  // Return the names of all fields
  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({
      "stress_xx","stress_yy","stress_zz","stress_yz","stress_zx","stress_xy",
      "strain_xx","strain_yy","strain_xy",
      "viscous_strain_xx","viscous_strain_yy","viscous_strain_xy"
    });
  }

  // Return the density
  Real density(void) { return rho; }
    
  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0;   // stress_xx
    state[1] = 0.0;   // stress_yy
    state[2] = 0.0;   // stress_zz
    state[3] = 0.0;   // stress_yz
    state[4] = 0.0;   // stress_zx
    state[5] = 0.0;   // stress_xy
    state[6] = 0.0;   // strain_xx
    state[7] = 0.0;   // strain_yy
    state[8] = 0.0;   // strain_xy
    state[9] = 0.0;   // viscous_strain_xx
    state[10] = 0.0;  // viscous_strain_yy
    state[11] = 0.0;  // viscous_strain_xy
  } // initialize()

  // Do we have to keep this?
  // Initialize variable material properties
  void initialize_variable_properties(Real (&x)[2], Real* state, double (*function_xy)(double,double)) {
    // Assign variable stiffness_scaling_factor as a function of initial spatial (x,y) coordinates
    // state[7] = function_xy(x[0],x[1]);
  } // initialize_variable_properties()

  // Update the material state using the current deformation gradient F
  void update(Real (&F)[2][2], Real &psi, Real* state, Real dt) {
    // Update at a material point.
    // Input: F (2x2), dt. Small strain: eps = sym(F - I).
    // Output: Updated states (and we potentially can have psi, elastic energy density.
    Real J = F[0][0]*F[1][1] - F[0][1]*F[1][0];
    if (J <= 0.0) {
      std::cout << "ERROR: negative Jacobian" << std::endl;
      exit(EXIT_FAILURE);
    }
    // 1) 2D small strain (engineering shear)
    Real strain_xx = F[0][0] - 1.0;
    Real strain_yy = F[1][1] - 1.0;
    Real strain_xy = 2*0.5*(F[0][1] + F[1][0]); // engineering gamma_xy                 
    Real strain[3] = { strain_xx, strain_yy, strain_xy };

    // 2) 2D deviatoric part 
    //    dev2(e) = e - 0.5*tr2(e)*I2, shears unchanged
    Real tr2 = strain_xx + strain_yy;
    Real dev_strain[3];
    dev_strain[0] = strain_xx - 0.5*tr2;
    dev_strain[1] = strain_yy - 0.5*tr2;
    dev_strain[2] = strain_xy;    
  } // update()

  // Conditionally load material history parameters from memory
  bool load_state(Real* state, Real* overflow_state)  { return false; }

  // Conditionally store material history parameters in memory
  bool store_state(Real* state, Real* overflow_state) { return false; }

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    field_data[0] = state[0]; // stress_xx
    field_data[1] = state[1]; // stress_yy
    field_data[2] = state[2]; // stress_zz
    field_data[3] = state[3]; // stress_yz
    field_data[4] = state[4]; // stress_zx
    field_data[5] = state[5]; // stress_xy
    field_data[6] = state[6]; // pressure
    field_data[7] = state[7]; // stiffness_scaling_factor
  }

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* ViscoElasticity */

#endif // VISCOELASTICITY_H
