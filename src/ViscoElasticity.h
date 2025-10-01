#ifndef VISCOELASTICITY_H
#define VISCOELASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"
#include "Dual.h"
#include "Fixed.h"
#include "Rational.h"
#include <limits>
#include "types.h"


template<class FixedE, class Ratio>
class ViscoElasticity {
 protected:
  Real rho;   // Mass density
  Real E;     // Young's modulus
  Real nu;    // Poisson's ratio
  Real mu;    // Shear modulus
  Real mu2;   // Twice the shear modulus
  Real lam;   // Lame parameter

  //Real kappa; // Bulk modulus
  //Real pmod;  // P-wave modulus
  
  // Viscous parameters
  Real mu_e;   // Shear modulus of the Maxwell element spring
  Real mu2_e;  // Twice the shear modulus of the Maxwell element spring
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

    // Get primary elastic constants ("equilibrium" stiffness parameters)
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


    //kappa = lam+0.5*mu2;
    //pmod  = lam+mu2;

    // Visco parameters
     // Get deviatoric/shear stiffness of Maxwell element spring (no bulk stiffness for Maxwell element)
    if (params.count("shear_modulus_Maxwell_element") > 0) {
      mu_e = params["shear_modulus_Maxwell_element"];
    } else {
      std::cout << "Missing material parameter: shear_modulus_Maxwell_element" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute shear constant for Maxwell element
    mu2_e    = 2*mu_e;

    // Gets either relatxation_time or viscosity
    // relatxation_time = viscosity / G (G = mu_e here)
    // If both provided, relaxation_time wins
    if (params.count("relaxation_time") > 0) {
      tau = params["relaxation_time"];
    } else {
      if (params.count("viscosity") > 0) {
        eta = params["viscosity"];
        tau = eta / mu_e;
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
 int num_state_vars(void) { return 15; }

  // Return the names of all fields
  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({
      "stress_xx","stress_yy","stress_zz","stress_yz","stress_zx","stress_xy",
      "strain_xx","strain_yy","strain_xy",
      "viscous_strain_xx","viscous_strain_yy","viscous_strain_xy",
      "viscous_strain_xx_dual","viscous_strain_yy_dual","viscous_strain_xy_dual"
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
    // store primal viscous strains as Fixed_E
    save_as_Real(FixedE(0.0), state[9]);   // viscous_strain_xx (primal)
    save_as_Real(FixedE(0.0), state[10]);  // viscous_strain_yy (primal)
    save_as_Real(FixedE(0.0), state[11]);  // viscous_strain_xy (primal)
    // store dual viscous strains as Fixed_E
    save_as_Real(FixedE(0.0), state[12]);  // viscous_strain_xx (dual)
    save_as_Real(FixedE(0.0), state[13]);  // viscous_strain_yy (dual)
    save_as_Real(FixedE(0.0), state[14]);  // viscous_strain_xy (dual)
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
 


  // 3) Viscous strain main update -- load dual/fixed stored values
  FixedE vs_xx_p, vs_xx_d;
  FixedE vs_yy_p, vs_yy_d;
  FixedE vs_xy_p, vs_xy_d;
  load_from_Real(state[9],  vs_xx_p); load_from_Real(state[12], vs_xx_d);
  load_from_Real(state[10], vs_yy_p); load_from_Real(state[13], vs_yy_d);
  load_from_Real(state[11], vs_xy_p); load_from_Real(state[14], vs_xy_d);

  Dual<FixedE> viscous_strain_xx(vs_xx_p, vs_xx_d);
  Dual<FixedE> viscous_strain_yy(vs_yy_p, vs_yy_d);
  Dual<FixedE> viscous_strain_xy(vs_xy_p, vs_xy_d);

  Real previous_strain_xx  =  state[6];
  Real previous_strain_yy  =  state[7];
  Real previous_strain_xy  =  state[8];
  // Define a vector version of previous strain for taking the dev part easier
  Real previous_strain[3] = { previous_strain_xx, previous_strain_yy,  previous_strain_xy }; // Previous step, strain
  Real dev_previous_strain[3];

  Real A = std::exp(-std::fabs(dt)/tau);
  // Use Rational for reversible operations
  Ratio A_rat(A);

    if (dt >= 0.0) {
      // First step of algorigthm
      viscous_strain_xx = viscous_strain_xx * A_rat;
      viscous_strain_yy = viscous_strain_yy * A_rat;
      viscous_strain_xy = viscous_strain_xy * A_rat;

      // Second step eps_prev = eps
      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      // Third step final calculation
      dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx(dev_previous_strain[0] * (1-A_rat), 0.0);
      Dual<FixedE> dyy(dev_previous_strain[1] * (1-A_rat), 0.0);
      Dual<FixedE> dxy(dev_previous_strain[2] * (1-A_rat), 0.0);
      viscous_strain_xx = viscous_strain_xx + dxx;
      viscous_strain_yy = viscous_strain_yy + dyy;
      viscous_strain_xy = viscous_strain_xy + dxy;

    } else {
      // First step of algorigthm
      dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx(dev_previous_strain[0] * (1-A_rat), 0.0);
      Dual<FixedE> dyy(dev_previous_strain[1] * (1-A_rat), 0.0);
      Dual<FixedE> dxy(dev_previous_strain[2] * (1-A_rat), 0.0);
           
      viscous_strain_xx = viscous_strain_xx - dxx;
      viscous_strain_yy = viscous_strain_yy - dyy;
      viscous_strain_xy = viscous_strain_xy - dxy;
     
      // Second step eps_prev = eps
      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      // Third step
      viscous_strain_xx = viscous_strain_xx / A_rat;
      viscous_strain_yy = viscous_strain_yy / A_rat;
      viscous_strain_xy = viscous_strain_xy / A_rat;
    } // end if dt

    // 4) Solve for stress with elastic part of strain
    // Equilibrium (long-time) stress: full isotropic Hooke in 2D
    Real stress_xx_eq = (lam + mu2)*strain_xx + lam*strain_yy;
    Real stress_yy_eq = lam* strain_xx + (lam + mu2)*strain_yy;
    Real stress_xy_eq = mu * strain_xy;

    // Maxwell branch (deviatoric only): sigma_M = 2*mu_e (deviatoric part of current strain - viscous part)
    // The reason why the previous_strain is used is that the current strain is poured into the previous_strain at this point. (after the update function)
    dev2(previous_strain,dev_previous_strain);
    Real dev_elastic_strain_xx = dev_previous_strain[0] - Real(viscous_strain_xx.first);
    Real dev_elastic_strain_yy = dev_previous_strain[1] - Real(viscous_strain_yy.first);
    Real dev_elastic_strain_xy = dev_previous_strain[2] - Real(viscous_strain_xy.first);

    Real stress_xx_Maxwell = mu2_e*(dev_elastic_strain_xx);
    Real stress_yy_Maxwell = mu2_e*(dev_elastic_strain_yy);
    Real stress_xy_Maxwell =  mu_e*(dev_elastic_strain_xy);

    // total stress = equilibrium + Maxwell deviatoric
    Real stress_xx = stress_xx_eq + stress_xx_Maxwell;
    Real stress_yy = stress_yy_eq + stress_yy_Maxwell;
    Real stress_xy = stress_xy_eq + stress_xy_Maxwell;

    // Becasue of being in 2D, explicitly set them to zero
    state[2]=0.0; state[3]=0.0; state[4]=0.0;

    // Update the states with final results
    state[0] = stress_xx;
    state[1] = stress_yy;
    state[5] = stress_xy;

    state[6] = previous_strain[0];
    state[7] = previous_strain[1];
    state[8] = previous_strain[2];

    save_as_Real(viscous_strain_xx.first,  state[9]);
    save_as_Real(viscous_strain_yy.first,  state[10]);
    save_as_Real(viscous_strain_xy.first,  state[11]);
    save_as_Real(viscous_strain_xx.second, state[12]);
    save_as_Real(viscous_strain_yy.second, state[13]);
    save_as_Real(viscous_strain_xy.second, state[14]);

    // elastic strain energy density
    // Because at this point we have the stress and total strain, we can compute the energy as 0.5*sigma:eps
    // Note that the shear strain is engineering shear strain
    // Elastic strain-energy density (springs only; dashpot excluded)
    Real psi_eq       = 0.5*(stress_xx_eq*strain_xx
                           + stress_yy_eq*strain_yy 
                           + stress_xy_eq*strain_xy);
    // elastic strain-energy density in Maxwell element
    Real psi_Maxwell  = 0.5*(stress_xx_Maxwell*dev_elastic_strain_xx
                           + stress_yy_Maxwell*dev_elastic_strain_yy
                           + stress_xy_Maxwell*dev_elastic_strain_xy);

    // Total elastic strain-energy density
    psi = psi_eq + psi_Maxwell;
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
    FixedE temp;
    load_from_Real(state[9] ,temp); field_data[9]  = Real(temp); // viscous_strain_xx
    load_from_Real(state[10],temp); field_data[10] = Real(temp); // viscous_strain_yy
    load_from_Real(state[11],temp); field_data[11] = Real(temp); // viscous_strain_xy
    load_from_Real(state[12],temp); field_data[12] = Real(temp); // viscous_strain_xx dual
    load_from_Real(state[13],temp); field_data[13] = Real(temp); // viscous_strain_yy dual
    load_from_Real(state[14],temp); field_data[14] = Real(temp); // viscous_strain_xy dual
  }

  // Return the initial sound speed
  //Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* ViscoElasticity */

#endif // VISCOELASTICITY_H
