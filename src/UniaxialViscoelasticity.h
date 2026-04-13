#ifndef UNIAXIAL_VISCOELASTICITY_H
#define UNIAXIAL_VISCOELASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include <array>
#include "Parameters.h"
#include "Dual.h"
#include "Fixed.h"
#include "PassPhase.h"
#include "Rational.h"
#include <limits>
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

  enum StateIndex {
    AXIAL_STRESS = 0,
    AXIAL_STRAIN = 1,
    VISCOUS_STRAIN = 2,
    DUAL_VISCOUS_STRAIN = 3,
    OVERFLOW_COUNTER = 4,
    DF_DTAU = 5,
    DF_DE = 6,
    LAMBDA_ADJOINT = 7,
    DUAL_LAMBDA_ADJOINT = 8
  };

 public:
  enum class ParameterId {
    RelaxationTime = 0,
    YoungsModulus = 1
  };

  struct ParameterGradient {
    std::array<Real,2> values;

    ParameterGradient() { values.fill(0.0); }

    Real& by_id(ParameterId id) { return values[int(id)]; }
    const Real& by_id(ParameterId id) const { return values[int(id)]; }
    Real& dTau(void) { return by_id(ParameterId::RelaxationTime); }
    Real& dE(void) { return by_id(ParameterId::YoungsModulus); }
    const Real& dTau(void) const { return by_id(ParameterId::RelaxationTime); }
    const Real& dE(void) const { return by_id(ParameterId::YoungsModulus); }
  };

  struct LocalForwardInput {
    Real eps_np1;
    Real q_n;
    Real dt;
  };

  struct LocalForwardOutput {
    Real A;
    Real sigma_np1;
    Real q_np1;
  };

  struct LocalBackwardInput {
    LocalForwardInput in;
    LocalForwardOutput out;
    Real bar_sigma_np1;
    Real bar_q_np1;
  };

  struct LocalBackwardOutput {
    Real bar_eps_np1;
    Real bar_q_n;
    ParameterGradient grad;
  };

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

  int num_state_vars(void) { return 9; }

  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({ "axial_stress",
                                      "axial_strain",
                                      "viscous_strain",
                                      "dual_viscous_strain",
                                      "overflow_counter",
                                      "df_dtau",
                                      "df_dE",
                                      "lambda_adjoint",
                                      "dual_lambda_adjoint" });
  }

  // Return the (mass per unit length) = (cross-sectional area) * (density)
  Real mass_per_unit_length(void) { return area*rho; }

  // Initialize the material state
  void initialize(Real* state) {
    state[AXIAL_STRESS] = 0.0;
    state[AXIAL_STRAIN] = 0.0;
    save_as_Real(FixedE(0.0), state[VISCOUS_STRAIN]);
    save_as_Real(FixedE(0.0), state[DUAL_VISCOUS_STRAIN]);
    state[OVERFLOW_COUNTER] = Real(0);
    state[DF_DTAU] = 0.0;
    state[DF_DE] = 0.0;
    save_as_Real(LambdaAdj(0.0), state[LAMBDA_ADJOINT]);
    save_as_Real(LambdaAdj(0.0), state[DUAL_LAMBDA_ADJOINT]);
  } // initialize()

  LocalForwardOutput forward_local(const LocalForwardInput& in) const {
    LocalForwardOutput out;
    out.A = std::exp(-std::fabs(in.dt)/tau);
    out.q_np1 = out.A*in.q_n + (1.0-out.A)*in.eps_np1;
    out.sigma_np1 = E*(in.eps_np1 - out.q_np1);
    return out;
  }

  LocalBackwardOutput backward_local(const LocalBackwardInput& in) const {
    LocalBackwardOutput out;
    const Real A = in.out.A;
    const Real de = in.in.eps_np1 - in.in.q_n;
    out.bar_eps_np1 = E*A*in.bar_sigma_np1 + (1.0-A)*in.bar_q_np1;
    out.bar_q_n = -E*A*in.bar_sigma_np1 + A*in.bar_q_np1;
    out.grad.dE() = A*de*in.bar_sigma_np1;
    out.grad.dTau() = de*(E*in.bar_sigma_np1 - in.bar_q_np1)*A*(std::fabs(in.in.dt)/(tau*tau));
    return out;
  }

  // Update the material state using the current stretch ratio
  void update(Real lambda, Real &psi, Real* state, Real dt, PassPhase phase) {
    // Check for element inversion
    if (lambda <= 0.0) {
      std::cout << "ERROR: non-positive stretch ratio" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Compute the axial (small) strain
    Real strain = lambda - 1.0;

    // Load the viscous strain and its dual from memory
    FixedE vs_p, vs_d;
    load_from_Real(state[VISCOUS_STRAIN], vs_p);
    load_from_Real(state[DUAL_VISCOUS_STRAIN], vs_d);
    Dual<FixedE> viscous_strain(vs_p, vs_d);

    // Accumulated gradient and adjoint lambda
    Real df_dtau_accum = state[DF_DTAU];
    Real df_dE_accum = state[DF_DE];
    LambdaAdj lambda_p, lambda_d;
    load_from_Real(state[LAMBDA_ADJOINT], lambda_p);
    load_from_Real(state[DUAL_LAMBDA_ADJOINT], lambda_d);
    Dual<LambdaAdj> lambda_adjoint(lambda_p, lambda_d);

    // Load the overflow counter from memory
    int overflow_counter = int(state[OVERFLOW_COUNTER]);

    // `state[AXIAL_STRAIN]` holds epsilon_n during forward and epsilon_{n+1} during reverse.
    Real previous_strain = state[AXIAL_STRAIN];

    const Real dt_abs = std::fabs(dt);
    const Real A_real = std::exp(-dt_abs/tau);
    Ratio A_rat(A_real);

    if (phase == PassPhase::Forward) {
      // Local forward kernel (explicit black-box view)
      LocalForwardInput fwd_in = { strain, Real(viscous_strain.first), dt_abs };
      LocalForwardOutput fwd_out = forward_local(fwd_in);
      (void)fwd_out;

      // Bit-reversible forward update of the history variable
      viscous_strain = viscous_strain * A_rat;
      previous_strain = strain;
      Dual<FixedE> dxx(previous_strain * (1-A_rat), 0.0);
      viscous_strain = viscous_strain + dxx;

      // Forward sweep terminal condition for reverse adjoint sweep.
      lambda_adjoint = Dual<LambdaAdj>(LambdaAdj(0.0), LambdaAdj(0.0));

      overflow_counter++;
    } else if (is_reverse_phase(phase)) {
      const bool run_adjoint = (phase == PassPhase::BackwardAdjoint);
      // For reverse mode, keep epsilon_{n+1} before reconstructing epsilon_n
      Real strain_n_plus_one = 0.0;
      if (run_adjoint) {
        strain_n_plus_one = previous_strain;
      }

      // Bit-reversible inverse update of the history variable
      Dual<FixedE> dxx(previous_strain * (1-A_rat), 0.0);
      viscous_strain = viscous_strain - dxx;
      previous_strain = strain;
      viscous_strain = viscous_strain / A_rat;

      if (run_adjoint) {
        // Adjoint state and gradient accumulation (objective seeding as in current visco examples)
        Real lambda_n_plus_one = Real(lambda_adjoint.first); // bar(q_{n+1})
        Real q_n = Real(viscous_strain.first);
        Real elastic_strain_n = strain - q_n;
        Real sigma_n = E * elastic_strain_n;

        // dL/dtau += bar(q_{n+1}) * dq_{n+1}/dtau
        ParameterGradient local_grad;
        local_grad.dTau() = lambda_n_plus_one * (q_n - strain_n_plus_one) * A_real * (dt_abs/(tau*tau));
        // d/dE of psi_n = 0.5*E*(eps_n-q_n)^2 with q treated as fixed local state input
        local_grad.dE() = 0.5*elastic_strain_n*elastic_strain_n;
        df_dtau_accum += local_grad.dTau();
        df_dE_accum += local_grad.dE();

        // Same style as viscous_strain reverse update: lambda_n from lambda_{n+1}
        // lambda_n = -sigma_n + A*lambda_{n+1}
        if constexpr (std::is_same<LambdaAdj, Real>::value) {
          lambda_adjoint = lambda_adjoint * A_real;
        } else {
          lambda_adjoint = lambda_adjoint * A_rat;
        }
        Dual<LambdaAdj> dlambda(LambdaAdj(-sigma_n), LambdaAdj(0.0));
        lambda_adjoint = lambda_adjoint + dlambda;
      }

      overflow_counter--;
    } else {
      std::cout << "UniaxialViscoelasticity::update received unsupported phase: "
                << pass_phase_name(phase) << std::endl;
      exit(EXIT_FAILURE);
    }

    state[AXIAL_STRAIN] = previous_strain;

    // Compute and store the axial stress
    Real elastic_strain = state[AXIAL_STRAIN] - Real(viscous_strain.first);
    Real stress = E * elastic_strain;
    state[AXIAL_STRESS] = stress;

    save_as_Real(viscous_strain.first,  state[VISCOUS_STRAIN]);
    save_as_Real(viscous_strain.second, state[DUAL_VISCOUS_STRAIN]);
    state[OVERFLOW_COUNTER] = Real(overflow_counter);
    state[DF_DTAU] = df_dtau_accum;
    state[DF_DE] = df_dE_accum;
    save_as_Real(lambda_adjoint.first,  state[LAMBDA_ADJOINT]);
    save_as_Real(lambda_adjoint.second, state[DUAL_LAMBDA_ADJOINT]);

    psi = 0.5 * stress * elastic_strain;
  }

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state) {
    if ((int(state[OVERFLOW_COUNTER]) == 0) && (overflow_state.size() >= 2)) {
      state[OVERFLOW_COUNTER] = Real(mat_overflow_limit);
      // Load in reverse order of storage
      state[DUAL_LAMBDA_ADJOINT] = overflow_state.back(); overflow_state.pop_back();
      state[DUAL_VISCOUS_STRAIN] = overflow_state.back(); overflow_state.pop_back();
    }
  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    if (int(state[OVERFLOW_COUNTER]) == mat_overflow_limit) {
      state[OVERFLOW_COUNTER] = Real(0);
      // Store in forward order to match reverse loading
      overflow_state.push_back(state[DUAL_VISCOUS_STRAIN]);
      overflow_state.push_back(state[DUAL_LAMBDA_ADJOINT]);
      state[DUAL_VISCOUS_STRAIN] = Real(0.0);
      state[DUAL_LAMBDA_ADJOINT] = Real(0.0);
    }
  }

  // Return the value of viscous_strain (Just for quick check, can be deleted later)
  // One could use the get_fields function
  Real get_state_variable(Real* state, std::string state_variable_name) {
    if (state_variable_name == "viscous_strain") {
      FixedE temp;
      load_from_Real(state[VISCOUS_STRAIN], temp);
      return Real(temp);
    } else if (state_variable_name == "df_dtau") {
      return state[DF_DTAU];
    } else if (state_variable_name == "df_dE") {
      return state[DF_DE];
    } else {
      return Real(0.0);
    }
  }

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    FixedE temp;
    LambdaAdj lambda_temp;
    field_data[0] = state[AXIAL_STRESS];
    field_data[1] = state[AXIAL_STRAIN];
    load_from_Real(state[VISCOUS_STRAIN], temp); field_data[2] = Real(temp);
    load_from_Real(state[DUAL_VISCOUS_STRAIN], temp); field_data[3] = Real(temp);
    field_data[4] = state[OVERFLOW_COUNTER];
    field_data[5] = state[DF_DTAU];
    field_data[6] = state[DF_DE];
    load_from_Real(state[LAMBDA_ADJOINT], lambda_temp); field_data[7] = Real(lambda_temp);
    load_from_Real(state[DUAL_LAMBDA_ADJOINT], lambda_temp); field_data[8] = Real(lambda_temp);
  }

  bool is_dead(Real*) { return false; }

}; /* UniaxialViscoelasticity */

#endif /* UNIAXIAL_VISCOELASTICITY_H */
