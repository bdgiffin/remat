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

  // Equlibrium spring parameters Standard Linear Solid (SLS) model
  Real E_eq;      // Equlibrium spring young's modulus
  Real nu_eq;     // Equlibrium spring Poisson's ratio
  Real lam_eq;    // Equlibrium spring lame parameter
  Real mu_eq;     // Equlibrium spring shear modulus
  Real mu2_eq;    // Twice the equlibrium spring shear modulus

  //Real kappa; // Bulk modulus
  //Real pmod;  // P-wave modulus
  // Viscous parameters
  Real tau;    // relaxation time
  Real eta;    // viscosity (if provided then tau = eta/mu would be computed)
 private:
  // Deviatoric projector in 2D
  inline void dev2(Real e[3], Real out[3]) {
    Real tr2 = e[0] + e[1];
    out[0] = e[0] - 0.5 * tr2; // eps_xx - tr/2
    out[1] = e[1] - 0.5 * tr2; // eps_yy - tr/2
    out[2] = e[2];             // gxy unchanged
  }

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
    
    // Get elastic constants for equilibrium spring
        if (params.count("youngs_modulus_equilibrium_spring") > 0) {
      E_eq = params["youngs_modulus_equilibrium_spring"];
    } else {
      std::cout << "Missing material parameter: youngs_modulus_equilibrium_spring" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (params.count("poissons_ratio_equilibrium_spring") > 0) {
      nu_eq = params["poissons_ratio_equilibrium_spring"];
    } else {
      std::cout << "Missing material parameter: poissons_ratio_equilibrium_spring" << std::endl;
      exit(EXIT_FAILURE);
    }
    

    // Compute derived elastic constants
    mu2   = E/(1.0+nu);
    mu    = 0.5*mu2;
    lam   = mu2*nu/(1.0-2.0*nu);

    // Compute derived elastic constants for equilibrium spring
    mu2_eq   = E_eq/(1.0+nu_eq);
    mu_eq    = 0.5*mu2_eq;
    lam_eq   = mu2_eq*nu_eq/(1.0-2.0*nu_eq);

    //kappa = lam+0.5*mu2;
    //pmod  = lam+mu2;

    // Visco parameters
    // Gets either relatxation_time or viscosity
    // relatxation_time = viscosity / G
    // If both provided, relaxation_time wins
    if (params.count("relaxation_time") > 0) {
      tau = params["relaxation_time"];
    } else {
      if (params.count("viscosity") > 0) {
        eta = params["viscosity"];
        tau = eta / mu;
      } else {      
        std::cout << "Missing viscous parameters: relaxation_time" << std::endl;
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
  void update(Real (&F)[2][2],Real &psi, Real* state, Real dt) {
    // Update at a material point.
    // Input: F (2x2), dt. Small strain: eps = sym(F - I).
    // Output: Updated states means updated viscous strain and stress values
    Real J = F[0][0]*F[1][1] - F[0][1]*F[1][0];
    if (J <= 0.0) {
      std::cout << "F = [[" << F[0][0] << ", " << F[0][1] << "], "<< "[" << F[1][0] << ", " << F[1][1] << "]]\n";
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
    //Real tr2 = strain_xx + strain_yy;
    //Real dev_strain[3];
    //dev_strain[0] = strain_xx - 0.5*tr2;
    //dev_strain[1] = strain_yy - 0.5*tr2;
    //dev_strain[2] = strain_xy;
    //dev2(strain,dev_strain);

    // 2) 2D deviatoric part of previous step strain 
    //    dev2(e) = e - 0.5*tr2(e)*I2, shears unchanged
    //Real tr2_previous = strain_previous[0] + strain_previous[1];
    Real dev_strain_previous[3];
    //dev_strain_previous[0] = strain_previous[0] - 0.5*tr2;
    //dev_strain_previous[1] = strain_previous[1] - 0.5*tr2;
    //dev_strain_previous[2] = strain_previous[2];
    //dev2(strain_previous,dev_strain_previous);


    // 3) Viscous strain main update
    Real viscous_strain[3]  = { state[9], state[10], state[11] }; // Previous step, viscous strain
    Real strain_previous[3] = { state[6], state[7],  state[8]  }; // Previous step, strain
    Real A = std::exp(-std::fabs(dt)/tau);

    if (dt >= 0.0) {
      // Solve for Viscous strain for step n+1
      //viscous_strain[0] = A*viscous_strain[0] + (1.0 - A)*dev_strain[0];
      //viscous_strain[1] = A*viscous_strain[1] + (1.0 - A)*dev_strain[1];
      //viscous_strain[2] = A*viscous_strain[2] + (1.0 - A)*dev_strain[2];

      // First step of algorigthm
      viscous_strain[0] *= A;
      viscous_strain[1] *= A;
      viscous_strain[2] *= A;

      // Second step eps_prev = eps
      strain_previous[0] = strain[0];
      strain_previous[1] = strain[1];
      strain_previous[2] = strain[2];

      // Third step final calculation
      dev2(strain_previous,dev_strain_previous);
      viscous_strain[0] += dev_strain_previous[0] * (1.0 - A);
      viscous_strain[1] += dev_strain_previous[1] * (1.0 - A);
      viscous_strain[2] += dev_strain_previous[2] * (1.0 - A);

    } else {
      // First step of algorigthm
      dev2(strain_previous,dev_strain_previous);
      viscous_strain[0] += - dev_strain_previous[0] * (1.0 - A);
      viscous_strain[1] += - dev_strain_previous[1] * (1.0 - A);
      viscous_strain[2] += - dev_strain_previous[2] * (1.0 - A);

      // Second step eps_prev = eps
      strain_previous[0] = strain[0];
      strain_previous[1] = strain[1];
      strain_previous[2] = strain[2];

      // Third step 
      viscous_strain[0] /= A;
      viscous_strain[1] /= A;
      viscous_strain[2] /= A;
    } // end if dt

    // 4) Solve for stress with elastic part of strain
    
    // Finding C for Maxwell 
    Real C11 = lam + mu2;
    Real C12 = lam;
    Real C66 = mu;
    // Finding C for equilibrium spring
    Real C11_eq = lam_eq + mu2_eq;
    Real C12_eq = lam_eq;
    Real C66_eq = mu_eq;

    Real elastic_strain_xx = strain[0] - viscous_strain[0];
    Real elastic_strain_yy = strain[1] - viscous_strain[1];
    Real elastic_strain_xy = strain[2] - viscous_strain[2];

    Real stress_xx = C11*elastic_strain_xx + C12*elastic_strain_yy + C11_eq*strain_xx + C12_eq*strain_yy;
    Real stress_yy = C12*elastic_strain_xx + C11*elastic_strain_yy + C12_eq*strain_xx + C11_eq*strain_yy;
    Real stress_xy = C66*elastic_strain_xy + C66_eq*strain_xy;
    Real stress[3] = { stress_xx, stress_yy, stress_xy };

    // Becasue of being in 2D, explicitly set them to zero
    state[2]=0.0; state[3]=0.0; state[4]=0.0;

    // Update the states with final results
    state[0] = stress[0];
    state[1] = stress[1];
    state[5] = stress[2];

    state[6] = strain[0];
    state[7] = strain[1];
    state[8] = strain[2];

    state[9]  = viscous_strain[0];
    state[10] = viscous_strain[1];
    state[11] = viscous_strain[2];
    psi = 1;
  } // update()

  // Conditionally load material history parameters from memory
  bool load_state(Real* state, Real* overflow_state)  { return false; }

  // Conditionally store material history parameters in memory
  bool store_state(Real* state, Real* overflow_state) { return false; }

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    field_data[0]  = state[0];  // stress_xx
    field_data[1]  = state[1];  // stress_yy
    field_data[2]  = state[2];  // stress_zz
    field_data[3]  = state[3];  // stress_yz
    field_data[4]  = state[4];  // stress_zx
    field_data[5]  = state[5];  // stress_xy
    field_data[6]  = state[6];  // strain_xx
    field_data[7]  = state[7];  // strain_yy
    field_data[8]  = state[8];  // strain_xy
    field_data[9]  = state[9];  // viscous_strain_xx
    field_data[10] = state[10]; // viscous_strain_yy
    field_data[11] = state[11]; // viscous_strain_xy
  }

  // Return the initial sound speed
  //Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* ViscoElasticity */

#endif // VISCOELASTICITY_H
