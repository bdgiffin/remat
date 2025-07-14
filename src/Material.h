#ifndef MATERIAL_H
#define MATERIAL_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"

class Material {
 protected:
  Real rho;   // Mass density
  Real E;     // Young's modulus
  Real nu;    // Poisson's ratio
  Real mu;    // Shear modulus
  Real mu2;   // Twice the shear modulus
  Real lam;   // Lame parameter
  Real kappa; // Bulk modulus
  Real pmod;  // P-wave modulus
 public:

  // Empty constructor
  Material(void) { }

  // Parameterized constructor
  Material(Parameters& params) {
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
  }
    
  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return 4; }

  // Return the density
  Real density(void) { return rho; }
    
  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0; // stress_xx
    state[1] = 0.0; // stress_yy
    state[2] = 0.0; // stress_xy
    state[3] = 0.0; // pressure
  } // initialize()
    
  // Update the material state using the current deformation gradient F
  void update(Real (&F)[2][2], Real &psi, Real* state, Real dt) {
    // Simple Neo-Hookean model:
    // J = det(F), Fbar = J^(-1/2) F, Bbar = Fbar*Fbar^T, Ibar = tr(Bbar), dev[Bbar] = Bbar - (Ibar/2)*1
    // psi = U(J) + W(Ibar)
    // tau = J*(dU/dJ)*1 + 2*(dW/dIbar)*dev[Bbar]
    // U = (kappa/2)*(J-1)^2 => dU/dJ    = (J-1)*kappa
    // W = (mu/2)*(Ibar - 2) => dW/dIbar = mu/2
    // sigma = (J-1)*kappa*1 + (mu/J)*dev[Bbar]

    // Compute the determinant of the deformation gradient J = det(F)
    Real J = F[0][0]*F[1][1] - F[0][1]*F[1][0];
    if (J <= 0.0) {
      std::cout << "ERROR: negative Jacobian" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute inverse of J
    Real invJ = 1.0/J;

    // Compute the left finger deformation tensor Bbar = Fbar*Fbar^T = (J^(-1/2))^2 (F*F^T)
    Real Bbar11 = invJ*(F[0][0]*F[0][0] + F[0][1]*F[0][1]);
    Real Bbar22 = invJ*(F[1][0]*F[1][0] + F[1][1]*F[1][1]);
    Real Bbar12 = invJ*(F[0][0]*F[1][0] + F[0][1]*F[1][1]);

    // Compute Ibar/2 = tr(Bbar)/2, and dev[Bbar] = Bbar - (Ibar/2)*1
    Real Ibar = Bbar11 + Bbar22;
    Bbar11 -= 0.5*Ibar;
    Bbar22 -= 0.5*Ibar;

    // Compute and store the Cauchy stress: sigma = (J-1)*kappa*1 + (mu/J)*dev[Bbar]
    Real p = kappa*(J-1.0);
    Real mu_bar = mu*invJ;
    state[0] = mu_bar*Bbar11 + p; // stress_xx
    state[1] = mu_bar*Bbar22 + p; // stress_yy
    state[2] = mu_bar*Bbar12;     // stress_xy
    state[3] = p;                 // pressure

    // Compute the strain energy density
    // Real U = 0.5*kappa*(J-1.0)*(J-1.0);
    // Real W = 0.5*mu*(Ibar - 2.0);
    psi = 0.5*p*(J-1.0) + 0.5*mu*(Ibar - 2.0);
    
  } // update()

  // Return the initial sound speed
  Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* Material */

#endif /* MATERIAL_H */
