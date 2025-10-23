#ifndef VISCOPLASTICITY_H
#define VISCOPLASTICITY_H

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
class ViscoPlasticity {
 protected:
  Real rho;   // Mass density
  Real E;     // Young's modulus
  Real nu;    // Poisson's ratio
  Real mu;    // Shear modulus
  Real mu2;   // Twice the shear modulus
  Real lam;   // Lame parameter
  Real yield_stress; // yield stress

  //Real kappa; // Bulk modulus
  //Real pmod;  // P-wave modulus

  // Viscoplastic parameters
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
  ViscoPlasticity(void) { }

  // Parameterized constructor
  ViscoPlasticity(Parameters& params) {
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

    // 2D relations
    //kappa = lam+0.5*mu2;
    //pmod  = lam+mu2;

    // Viscoplastic parameters
     // Get deviatoric/shear stiffness of Maxwell element spring (no bulk stiffness for Maxwell element)
    if (params.count("shear_modulus_Maxwell_element") > 0) {
      mu_e = params["shear_modulus_Maxwell_element"];
    } else {
      std::cout << "Missing material parameter: shear_modulus_Maxwell_element" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (params.count("yield_stress") > 0) {
      yield_stress = params["yield_stress"];
    } else {
      std::cout << "Missing material parameter: yield_stress" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute shear constant for Maxwell element
    mu2_e    = 2.0*mu_e;

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
        std::cout << "Missing viscoplastic parameters: relaxation_time" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (tau<=0.0) {
      std::cout << "ViscoPlasticity: need a positive relaxation_time" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
  // Return the number of state variables for allocation purposes
 int num_state_vars(void) { return 16; }

  // Return the names of all fields
  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({
      "stress_xx","stress_yy","stress_zz","stress_yz","stress_zx","stress_xy",
      "strain_xx","strain_yy","strain_xy",
      "plastic_strain_xx","plastic_strain_yy","plastic_strain_xy",
      "plastic_strain_xx_dual","plastic_strain_yy_dual","plastic_strain_xy_dual",
      "stiffness_scaling_factor"
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
    // store primal plastic strains as Fixed_E
    save_as_Real(FixedE(0.0), state[9]);   // plastic_strain_xx (primal)
    save_as_Real(FixedE(0.0), state[10]);  // plastic_strain_yy (primal)
    save_as_Real(FixedE(0.0), state[11]);  // plastic_strain_xy (primal)
    // store dual plastic strains as Fixed_E
    save_as_Real(FixedE(0.0), state[12]);  // plastic_strain_xx (dual)
    save_as_Real(FixedE(0.0), state[13]);  // plastic_strain_yy (dual)
    save_as_Real(FixedE(0.0), state[14]);  // plastic_strain_xy (dual)
    state[15] = 1.0; // stiffness_scaling_factor
  } // initialize()

  // Do we have to keep this?
  // Initialize variable material properties
  void initialize_variable_properties(Real (&x)[2], Real* state, double (*function_xy)(double,double)) {
    // Assign variable stiffness_scaling_factor as a function of initial spatial (x,y) coordinates
    state[15] = function_xy(x[0],x[1]);
  } // initialize_variable_properties()

  // Update the material state using the current deformation gradient F
  void update(Real (&F)[2][2],Real &psi, Real* state, Real dt) {
    // Update at a material point.
    // Input: F (2x2), dt. Small strain: eps = sym(F - I).
    // Output: Updated states means updated plastic strain and stress values
    Real J = F[0][0]*F[1][1] - F[0][1]*F[1][0];
    if (J <= 0.0) {
      std::cout << "F = [[" << F[0][0] << ", " << F[0][1] << "], "<< "[" << F[1][0] << ", " << F[1][1] << "]]\n";
      std::cout << "ERROR: negative Jacobian" << std::endl;
      exit(EXIT_FAILURE);
    }

    // 1) 2D small strain (engineering shear)
    Real strain_xx = F[0][0] - 1.0;
    Real strain_yy = F[1][1] - 1.0;
    Real strain_xy = 2.0*0.5*(F[0][1] + F[1][0]); // engineering gamma_xy
             
 


  // 3) plastic strain main update -- load dual/fixed stored values
  FixedE ps_xx_p, ps_xx_d;
  FixedE ps_yy_p, ps_yy_d;
  FixedE ps_xy_p, ps_xy_d;
  load_from_Real(state[9],  ps_xx_p); load_from_Real(state[12], ps_xx_d);
  load_from_Real(state[10], ps_yy_p); load_from_Real(state[13], ps_yy_d);
  load_from_Real(state[11], ps_xy_p); load_from_Real(state[14], ps_xy_d);

  Dual<FixedE> plastic_strain_xx(ps_xx_p, ps_xx_d);
  Dual<FixedE> plastic_strain_yy(ps_yy_p, ps_yy_d);
  Dual<FixedE> plastic_strain_xy(ps_xy_p, ps_xy_d);

  Real previous_strain_xx  =  state[6];
  Real previous_strain_yy  =  state[7];
  Real previous_strain_xy  =  state[8];
  // Define a vector version of previous strain for taking the dev part easier
  Real previous_strain[3] = { previous_strain_xx, previous_strain_yy,  previous_strain_xy }; // Previous step, strain
  //Real dev_previous_strain[3];
  Real sigma_trial_xx =    E * (strain_xx - Real(plastic_strain_xx.first));
  Real f_y_trial = std::fabs(sigma_trial_xx) - (yield_stress);
  //std::cout << "sigma_trial_xx = " << sigma_trial_xx << std::endl;
  //std::cout << "f_y_trial = " << f_y_trial << std::endl;
  if (f_y_trial <= 0) {
  Real A = std::exp(-std::fabs(dt)/tau);
  // Use Rational for reversible operations
  Ratio A_rat(A);
  Real n = std::copysign(1.0,sigma_trial_xx); 
    if (dt >= 0.0) {
      // First step of algorigthm
      plastic_strain_xx = plastic_strain_xx * A_rat;
      plastic_strain_yy = plastic_strain_yy * A_rat;
      plastic_strain_xy = plastic_strain_xy * A_rat;

      // Second step eps_prev = eps
      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      // Third step final calculation
      //dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx((previous_strain[0]-n*sigma_trial_xx/E) * (1-A_rat), 0.0);
      //Dual<FixedE> dyy(previous_strain[1] * (1-A_rat), 0.0);
      //Dual<FixedE> dxy(previous_strain[2] * (1-A_rat), 0.0);
      plastic_strain_xx = plastic_strain_xx + dxx;
      //plastic_strain_yy = plastic_strain_yy + dyy;
      //plastic_strain_xy = plastic_strain_xy + dxy;

    } else {
      // First step of algorigthm
      //dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx((previous_strain[0]-n*sigma_trial_xx/E) * (1-A_rat), 0.0);
      //Dual<FixedE> dyy(previous_strain[1] * (1-A_rat), 0.0);
      //Dual<FixedE> dxy(previous_strain[2] * (1-A_rat), 0.0);

      plastic_strain_xx = plastic_strain_xx - dxx;
      //plastic_strain_yy = plastic_strain_yy - dyy;
      //plastic_strain_xy = plastic_strain_xy - dxy;

      // Second step eps_prev = eps
      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      // Third step
      plastic_strain_xx = plastic_strain_xx / A_rat;
      plastic_strain_yy = plastic_strain_yy / A_rat;
      plastic_strain_xy = plastic_strain_xy / A_rat;
    } // end if dt
  } else {
    previous_strain[0] = strain_xx;
    previous_strain[1] = strain_yy;
    previous_strain[2] = strain_xy;
    
  }

    // 4) Solve for stress with elastic part of strain
    Real stiffness_scaling_factor = state[15];

    // Scaled equilibrium stiffnesses
    Real mu2_scaled = stiffness_scaling_factor * mu2;
    Real mu_scaled  = stiffness_scaling_factor * mu;
    // Scaled first Lame parameter (Watch out this is only exact for 2D)
    Real lam_scaled = stiffness_scaling_factor * lam;
    // Apply scaling to the Maxwell shear stiffness as well
    Real mu2_e_scaled = stiffness_scaling_factor * mu2_e;
    Real mu_e_scaled  = stiffness_scaling_factor * mu_e;

    // Equilibrium (long-time) stress: full isotropic Hooke in 2D
    Real stress_xx_eq = (lam_scaled + mu2_scaled)*strain_xx + lam_scaled*strain_yy;
    Real stress_yy_eq = lam_scaled* strain_xx + (lam_scaled + mu2_scaled)*strain_yy;
    Real stress_xy_eq = mu_scaled * strain_xy;

    // Maxwell branch (deviatoric only): sigma_M = 2*mu_e (deviatoric part of current strain - plastic part)
    // The reason why the previous_strain is used is that the current strain is poured into the previous_strain at this point. (after the update function)
    //dev2(previous_strain,dev_previous_strain);
    Real elastic_strain_xx = previous_strain[0] - Real(plastic_strain_xx.first);
    Real elastic_strain_yy = previous_strain[1] - Real(plastic_strain_yy.first);
    Real elastic_strain_xy = previous_strain[2] - Real(plastic_strain_xy.first);

    Real stress_xx_Maxwell = mu2_e_scaled*(elastic_strain_xx);
    Real stress_yy_Maxwell = mu2_e_scaled*(elastic_strain_yy);
    Real stress_xy_Maxwell = mu_e_scaled*(elastic_strain_xy);

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

    save_as_Real(plastic_strain_xx.first,  state[9]);
    save_as_Real(plastic_strain_yy.first,  state[10]);
    save_as_Real(plastic_strain_xy.first,  state[11]);
    save_as_Real(plastic_strain_xx.second, state[12]);
    save_as_Real(plastic_strain_yy.second, state[13]);
    save_as_Real(plastic_strain_xy.second, state[14]);

    // elastic strain energy density
    // Because at this point we have the stress and total strain, we can compute the energy as 0.5*sigma:eps
    // Note that the shear strain is engineering shear strain
    // Elastic strain-energy density (springs only; dashpot excluded)
    Real psi_eq       = 0.5*(stress_xx_eq*strain_xx
                           + stress_yy_eq*strain_yy 
                           + stress_xy_eq*strain_xy);
    // elastic strain-energy density in Maxwell element
    Real psi_Maxwell  = 0.5*(stress_xx_Maxwell*elastic_strain_xx
                           + stress_yy_Maxwell*elastic_strain_yy
                           + stress_xy_Maxwell*elastic_strain_xy);

    // Total elastic strain-energy density
    psi = psi_eq + psi_Maxwell;
  } // update()

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state) { (void)overflow_state; }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) { (void)overflow_state; }

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
    load_from_Real(state[9] ,temp); field_data[9]  = Real(temp); // plastic_strain_xx
    load_from_Real(state[10],temp); field_data[10] = Real(temp); // plastic_strain_yy
    load_from_Real(state[11],temp); field_data[11] = Real(temp); // plastic_strain_xy
    load_from_Real(state[12],temp); field_data[12] = Real(temp); // plastic_strain_xx dual
    load_from_Real(state[13],temp); field_data[13] = Real(temp); // plastic_strain_yy dual
    load_from_Real(state[14],temp); field_data[14] = Real(temp); // plastic_strain_xy dual
    field_data[15] = state[15]; // stiffness_scaling_factor
  }

  // Return the initial sound speed
  //Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* ViscoPlasticity */

#endif // VISCOPLASTICITY_H
