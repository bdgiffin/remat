#ifndef UNIAXIAL_VISCOELASTICITY_H
#define UNIAXIAL_VISCOELASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"
#include "ConstitutiveAdjoint.h"
#include "Dual.h"
#include "Fixed.h"
#include "Rational.h"
#include <limits>
#include <string>
#include <type_traits>
#include "types.h"


template<class FixedE, class Ratio, class LambdaAdj = FixedE>
class UniaxialViscoelasticity {
 protected:
  Real rho;   // Mass density
  Real area;  // Cross-sectional area
  Real E;     // Young's modulus
  Real eta;   // Viscosity
  Real tau;   // Relaxation time
  int  mat_overflow_limit = std::numeric_limits<int>::max();

 public:
  UniaxialViscoelasticity(void) { }

  UniaxialViscoelasticity(Parameters& params) {
    // Get material density
    if  (params.count("truss_density") > 0) {
      rho = params["truss_density"];
    } else if (params.count("density") > 0) {
      rho = params["density"];
    } else {
      std::cout << "Missing uniaxial viscoelasticity parameter: density / truss_density" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get cross-sectional area
    if (params.count("area") > 0) {
      area = params["area"];
    } else {
      std::cout << "Missing uniaxial viscoelastic parameter: area" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get primary elastic constants
    if  (params.count("truss_youngs_modulus") > 0) {
      E = params["truss_youngs_modulus"];
    } else if (params.count("youngs_modulus") > 0) {
      E = params["youngs_modulus"];
    } else {
      std::cout << "Missing uniaxial viscoelasticity parameter: youngs_modulus / truss_youngs_modulus" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Verify that the Young's modulus is a strictly positive value:
    if (E <= 0.0) {
      std::cout << "Invalid uniaxial viscoplasticity parameter specified: (youngs_modulus / truss_youngs_modulus) must be a strictly positive value" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get viscosity or relaxation_time
    if (params.count("relaxation_time") > 0) {
      tau = params["relaxation_time"];
      if (tau <= 0.0) {
        std::cout << "Uniaxial viscoelasticity requires a positive relaxation_time" << std::endl;
        exit(EXIT_FAILURE);
      }
      eta = tau * E;
    } else if (params.count("viscosity") > 0) {
      eta = params["viscosity"];
      if (eta <= 0.0) {
        std::cout << "Uniaxial viscoelasticity requires a positive viscosity" << std::endl;
        exit(EXIT_FAILURE);
      }
      tau = eta / E;
    } else {
      std::cout << "Missing uniaxial viscoelastic parameter: relaxation_time or viscosity" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get overflow limit for material history variables
    if (params.count("mat_overflow_limit") > 0) { mat_overflow_limit = int(params["mat_overflow_limit"]); }
  }

  int num_state_vars(void) { return 7; }

  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({ "axial_stress",
                                      "axial_strain",
                                      "viscous_strain",
                                      "dual_viscous_strain",
                                      "overflow_counter",
                                      "dL_dparam_relaxation_time",
                                      "lambda_adjoint" });
  }

  // Return the (mass per unit length) = (cross-sectional area) * (density)
  Real mass_per_unit_length(void) { return area*rho; }

  // Return Young's modulus used in the constitutive update
  Real youngs_modulus(void) { return E; }

  Real axial_stress(const Real* state) const { return state[0]; }

  Real axial_tangent_modulus(const Real* /*state*/) const { return E; }

  void add_reverse_stress_seed(Real* state, Real dL_dstress) {
    LambdaAdj lambda_p;
    load_from_Real(state[6], lambda_p);
    lambda_p = lambda_p + LambdaAdj(-E*dL_dstress);
    save_as_Real(lambda_p, state[6]);
  }

  // Generic parameter metadata for black-box gradient aggregation in System.
  int num_params(void) const { return 1; }

  const char* param_name(int param_id) const {
    if (param_id == 0) { return "relaxation_time"; }
    return nullptr;
  }

  Real param_value(int param_id) const {
    if (param_id == 0) { return tau; }
    return 0.0;
  }

  Real get_param_gradient(const Real* state, int param_id) const {
    if (param_id == 0) { return state[5]; }
    return 0.0;
  }

  // Adapter methods for the generic constitutive interface.
  void forward(const MaterialPointKinematics<Real>& kin, Real* state, MaterialForwardOutput<Real>& out) {
    update(kin.primary_measure,out.psi,state,kin.dt);
    out.stress = state[0];
  }

  void reverse(const MaterialPointKinematics<Real>& kin, Real* state, const MaterialReverseSeed<Real>& seed, MaterialReverseOutput<Real>& out) {
    out.dL_dmeasure = 0.0;
    out.dL_dparams.assign(num_params(),0.0);
    // Assumption: reverse path remains stateful in update(dt<0) for this model version.
    add_reverse_stress_seed(state,seed.dL_dstress);
    Real psi = 0.0;
    update(kin.primary_measure,psi,state,kin.dt);
    out.stress = state[0];
    out.psi = psi;
    out.dL_dparams[0] = state[5];
  }

  // Initialize the material state
  void initialize(Real* state) {
    state[0] = 0.0; // axial_stress
    state[1] = 0.0; // axial_strain
    save_as_Real(FixedE(0.0), state[2]); // viscous_strain
    save_as_Real(FixedE(0.0), state[3]); // dual_viscous_strain
    state[4] = Real(0); // overflow_counter
    state[5] = 0.0; // accumulated sensitivity dL/d(relaxation_time)
    save_as_Real(LambdaAdj(0.0), state[6]); // lambda_adjoint
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

    // Load the viscous strain and its dual from memory
    FixedE vs_p, vs_d;
    load_from_Real(state[2], vs_p);
    load_from_Real(state[3], vs_d);
    Dual<FixedE> viscous_strain(vs_p, vs_d);

    // Accumulated gradient and adjoint lambda
    Real dL_drelaxation_time_accum = state[5];
    LambdaAdj lambda_p;
    load_from_Real(state[6], lambda_p);
    LambdaAdj lambda_adjoint = lambda_p;

    // Load the overflow counter from memory
    int overflow_counter = int(state[4]);

    // Conditionally update the strain history variable
    // if (dt > 0.0) { state[1] = strain; } // strain_prev <- strain
    // Dual<FixedE> total_strain(FixedE(state[1]), FixedE(0.0));
    Real previous_strain = state[1];

    Real dt_abs = std::fabs(dt);
    // A = exp(-|dt|/tau)
    Real A_real = std::exp(-dt_abs/tau);
    Ratio A_rat(A_real);

    if (dt >= 0.0) {
      viscous_strain = viscous_strain * A_rat;
      previous_strain = strain;
      Dual<FixedE> dxx (previous_strain * (1-A_rat), 0.0);
      viscous_strain = viscous_strain + dxx;

      // Forward sweep terminal condition for reverse adjoint sweep.
      lambda_adjoint = LambdaAdj(0.0);

      overflow_counter++;
    } else {
      // For adjoint mode, we need to keep the strain value of n+1 step
      Real strain_n_plus_one = previous_strain;

      // Simple viscous strain update
      Dual<FixedE> dxx (previous_strain * (1-A_rat), 0.0);
      viscous_strain = viscous_strain - dxx;
      previous_strain = strain;
      viscous_strain = viscous_strain / A_rat;

      // Adjoint variable update (lambda) and gradient assembly
      // lambda_n = -sigma_n + A_n * lambda_{n+1}
      // dL/dtau += lambda_{n+1} * (viscous_strain_n - eps_{n+1}) * A_n * dt / tau^2
      // One assumption here is that the Psi function is not dependent on the last state (No final terminal in the loss function f)
      // So that we can initialize at the end lambda_{N} = 0.0 and this happens in the intialize() function above.
      Real lambda_n_plus_one = Real(lambda_adjoint); // lambda_{n+1} stored in state before step
      
      Real vs_prev = Real(viscous_strain.first); 
      Real elastic_strain_prev = strain - vs_prev; 
      Real sigma_n = E * elastic_strain_prev;

      dL_drelaxation_time_accum += lambda_n_plus_one * (vs_prev - strain_n_plus_one) * A_real * (dt_abs/(tau*tau));

      // Same style as viscous_strain reverse update: lambda_n from lambda_{n+1}
      // lambda_n = -sigma_n + A*lambda_{n+1}
      if constexpr (std::is_same<LambdaAdj, Real>::value) {
        lambda_adjoint = lambda_adjoint * A_real;
      } else {
        lambda_adjoint = lambda_adjoint * A_rat;
      }
      lambda_adjoint = lambda_adjoint + LambdaAdj(-sigma_n);

      overflow_counter--;
    }

    state[1] = previous_strain;

    // Compute and store the axial stress
    Real elastic_strain = state[1] - Real(viscous_strain.first);
    Real stress = E * elastic_strain;
    state[0] = stress;

    save_as_Real(viscous_strain.first,  state[2]);
    save_as_Real(viscous_strain.second, state[3]);
    state[4] = Real(overflow_counter);
    state[5] = dL_drelaxation_time_accum;
    save_as_Real(lambda_adjoint, state[6]);

    psi = 0.5 * stress * elastic_strain;
  }

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state) {
    if ((int(state[4]) == 0) && (overflow_state.size() >= 1)) {
      state[4] = Real(mat_overflow_limit);
      state[3] = overflow_state.back(); overflow_state.pop_back(); // dual_viscous_strain
    }
  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    if (int(state[4]) == mat_overflow_limit) {
      state[4] = Real(0);
      overflow_state.push_back(state[3]); // dual_viscous_strain
      state[3] = Real(0.0);
    }
  }

  // Return the value of viscous_strain (Just for quick check, can be deleted later)
  // One could use the get_fields function
  Real get_state_variable(Real* state, std::string state_variable_name) {
    if (state_variable_name == "viscous_strain") {
      FixedE temp;
      load_from_Real(state[2], temp);
      return Real(temp);
    } else {
      return Real(0.0);
    }
  }

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    FixedE temp;
    LambdaAdj lambda_temp;
    field_data[0] = state[0]; // axial_stress
    field_data[1] = state[1]; // axial_strain
    load_from_Real(state[2], temp); field_data[2] = Real(temp); // viscous_strain
    load_from_Real(state[3], temp); field_data[3] = Real(temp); // dual_viscous_strain
    field_data[4] = state[4]; // overflow_counter
    field_data[5] = state[5]; // dL_dparam_relaxation_time accumulation
    load_from_Real(state[6], lambda_temp); field_data[6] = Real(lambda_temp); // lambda_adjoint
  }

  bool is_dead(Real*) { return false; }

}; /* UniaxialViscoelasticity */

#endif /* UNIAXIAL_VISCOELASTICITY_H */
