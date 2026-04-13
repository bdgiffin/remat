#ifndef UNIAXIAL_VISCOELASTICITY_H
#define UNIAXIAL_VISCOELASTICITY_H

#include <math.h>
#include <cmath>
#include <iostream>
#include <stdlib.h> // exit
#include <limits>
#include <type_traits>
#include <vector>

#include "Parameters.h"
#include "Dual.h"
#include "Fixed.h"
#include "Rational.h"
#include "types.h"
#include "AdjointFramework.h"

template<class FixedE, class Ratio, class LambdaAdj = FixedE>
class UniaxialViscoelasticity {
 public:
  enum class MaterialUpdateMode {
    Forward,
    RematBackward,
    AdjointBackward
  };

  struct MaterialPrimalState {
    Real axial_stress = 0.0;
    Real axial_strain = 0.0; // Stores epsilon_n for the current state n
    Dual<FixedE> viscous_strain = Dual<FixedE>(FixedE(0.0),FixedE(0.0));
    int overflow_counter = 0;
  };

  struct MaterialAdjointState {
    // Local adjoint multiplier associated with viscous history variable.
    Dual<LambdaAdj> lambda_history = Dual<LambdaAdj>(LambdaAdj(0.0),LambdaAdj(0.0));
  };

  struct GradientAccumulator {
    Real df_dtau = 0.0;
    Real df_dE = 0.0;
  };

  struct MaterialStepInput {
    Real strain_n = 0.0;    // epsilon_n
    Real strain_np1 = 0.0;  // epsilon_{n+1}
    Real dt = 0.0;
    MaterialUpdateMode mode = MaterialUpdateMode::Forward;
  };

  struct LocalAdjointSeed {
    // Seed for stress sensitivity at step n.
    Real bar_sigma_n = 0.0;

    // Optional seed for stress sensitivity at step n+1.
    Real bar_sigma_np1 = 0.0;

    // Optional direct contribution for explicit objective dependence on E.
    Real direct_dE = 0.0;
  };

  struct MaterialForwardResponse {
    Real stress = 0.0;
    Real energy = 0.0;
  };

  struct MaterialAdjointResponse {
    Real bar_strain_n = 0.0;
    Real bar_strain_np1 = 0.0;
    Real grad_dtau_increment = 0.0;
    Real grad_dE_increment = 0.0;
  };

  struct MaterialStepRecord {
    Real dt = 0.0;
    Real strain_n = 0.0;
    Real strain_np1 = 0.0;
  };

  class StepTape {
   public:
    void clear() { records.clear(); }
    void push(const MaterialStepRecord& rec) { records.push_back(rec); }
    size_t size() const { return records.size(); }
    const MaterialStepRecord& at(size_t i) const { return records[i]; }

   private:
    std::vector<MaterialStepRecord> records;
  };

 protected:
  Real rho;   // Mass density
  Real area;  // Cross-sectional area
  Real E;     // Young's modulus
  Real eta;   // Viscosity
  Real tau;   // Relaxation time
  int  mat_overflow_limit = std::numeric_limits<int>::max();
  bool m_enable_adjoint = true;
  bool m_has_external_seed = false;
  LocalAdjointSeed m_external_seed;

  enum StateIndex {
    kAxialStress = 0,
    kAxialStrain = 1,
    kViscousStrain = 2,
    kDualViscousStrain = 3,
    kOverflowCounter = 4,
    kDfDtau = 5,
    kDfDE = 6,
    kLambdaAdjoint = 7,
    kDualLambdaAdjoint = 8,
    kStateSize = 9
  };

  MaterialPrimalState load_primal_state(const Real* state) const {
    MaterialPrimalState primal;
    primal.axial_stress = state[kAxialStress];
    primal.axial_strain = state[kAxialStrain];

    FixedE vs_primal, vs_dual;
    load_from_Real(state[kViscousStrain],vs_primal);
    load_from_Real(state[kDualViscousStrain],vs_dual);
    primal.viscous_strain = Dual<FixedE>(vs_primal,vs_dual);

    primal.overflow_counter = int(state[kOverflowCounter]);
    return primal;
  }

  void store_primal_state(const MaterialPrimalState& primal, Real* state) const {
    state[kAxialStress] = primal.axial_stress;
    state[kAxialStrain] = primal.axial_strain;
    save_as_Real(primal.viscous_strain.first, state[kViscousStrain]);
    save_as_Real(primal.viscous_strain.second,state[kDualViscousStrain]);
    state[kOverflowCounter] = Real(primal.overflow_counter);
  }

  MaterialAdjointState load_adjoint_state(const Real* state) const {
    MaterialAdjointState adj;
    LambdaAdj lam_primal, lam_dual;
    load_from_Real(state[kLambdaAdjoint],lam_primal);
    load_from_Real(state[kDualLambdaAdjoint],lam_dual);
    adj.lambda_history = Dual<LambdaAdj>(lam_primal,lam_dual);
    return adj;
  }

  void store_adjoint_state(const MaterialAdjointState& adj, Real* state) const {
    save_as_Real(adj.lambda_history.first,state[kLambdaAdjoint]);
    save_as_Real(adj.lambda_history.second,state[kDualLambdaAdjoint]);
  }

  GradientAccumulator load_gradients(const Real* state) const {
    GradientAccumulator grads;
    grads.df_dtau = state[kDfDtau];
    grads.df_dE = state[kDfDE];
    return grads;
  }

  void store_gradients(const GradientAccumulator& grads, Real* state) const {
    state[kDfDtau] = grads.df_dtau;
    state[kDfDE] = grads.df_dE;
  }

  MaterialForwardResponse make_response(const MaterialPrimalState& primal) const {
    MaterialForwardResponse response;
    const Real elastic_strain = primal.axial_strain - Real(primal.viscous_strain.first);
    response.stress = primal.axial_stress;
    response.energy = 0.5*primal.axial_stress*elastic_strain;
    return response;
  }

  void check_stretch(Real lambda) const {
    if (lambda <= 0.0) {
      std::cout << "ERROR: non-positive stretch ratio" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  LocalAdjointSeed consume_objective_seed(const MaterialPrimalState& primal_n) {
    if (m_has_external_seed) {
      m_has_external_seed = false;
      return m_external_seed;
    }
    return make_objective_seed(primal_n);
  }

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

    // Verify that the Young's modulus is a strictly positive value
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
      eta = tau*E;
    } else if (params.count("viscosity") > 0) {
      eta = params["viscosity"];
      if (eta <= 0.0) {
        std::cout << "Uniaxial viscoelasticity requires a positive viscosity" << std::endl;
        exit(EXIT_FAILURE);
      }
      tau = eta/E;
    } else {
      std::cout << "Missing uniaxial viscoelastic parameter: relaxation_time or viscosity" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Get overflow limit for material history variables
    if (params.count("mat_overflow_limit") > 0) {
      mat_overflow_limit = int(params["mat_overflow_limit"]);
    }
  }

  int num_state_vars(void) { return kStateSize; }

  // Accessors and control hooks used by the generalized driver.
  Real youngs_modulus(void) const { return E; }
  int adjoint_support_level(void) { return int(AdjointSupportLevel::Supported); }
  const char* adjoint_support_status(void) { return "supported"; }
  void set_adjoint_enabled(bool enabled) { m_enable_adjoint = enabled; }
  void set_external_objective_seed(const LocalAdjointSeed& seed) {
    m_external_seed = seed;
    m_has_external_seed = true;
  }
  void clear_external_objective_seed(void) { m_has_external_seed = false; }
  void reset_adjoint_state(Real* state) {
    state[kDfDtau] = 0.0;
    state[kDfDE] = 0.0;
    save_as_Real(LambdaAdj(0.0),state[kLambdaAdjoint]);
    save_as_Real(LambdaAdj(0.0),state[kDualLambdaAdjoint]);
  }

  MaterialPrimalState read_primal_state(const Real* state) const { return load_primal_state(state); }
  MaterialAdjointState read_adjoint_state(const Real* state) const { return load_adjoint_state(state); }
  GradientAccumulator read_gradient_state(const Real* state) const { return load_gradients(state); }
  void write_primal_state(const MaterialPrimalState& primal, Real* state) const { store_primal_state(primal,state); }
  void write_adjoint_state(const MaterialAdjointState& adjoint, Real* state) const { store_adjoint_state(adjoint,state); }
  void write_gradient_state(const GradientAccumulator& gradients, Real* state) const { store_gradients(gradients,state); }

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
    state[kAxialStress] = 0.0;
    state[kAxialStrain] = 0.0;
    save_as_Real(FixedE(0.0),state[kViscousStrain]);
    save_as_Real(FixedE(0.0),state[kDualViscousStrain]);
    state[kOverflowCounter] = Real(0);
    state[kDfDtau] = 0.0;
    state[kDfDE] = 0.0;
    save_as_Real(LambdaAdj(0.0),state[kLambdaAdjoint]);
    save_as_Real(LambdaAdj(0.0),state[kDualLambdaAdjoint]);
  }

  MaterialForwardResponse forward_update(const MaterialStepInput& input,
                                         MaterialPrimalState& primal,
                                         MaterialAdjointState& adjoint) const {
    const Real dt_abs = std::fabs(input.dt);
    const Real A_real = std::exp(-dt_abs/tau);
    Ratio A_rat(A_real);

    // ev_{n+1} = A*ev_n + (1-A)*eps_{n+1}
    Dual<FixedE> dxx(input.strain_np1*(1-A_rat),0.0);
    primal.viscous_strain = primal.viscous_strain*A_rat;
    primal.viscous_strain = primal.viscous_strain + dxx;

    primal.axial_strain = input.strain_np1;

    // sigma_{n+1} = E*(eps_{n+1} - ev_{n+1})
    const Real elastic_strain = primal.axial_strain - Real(primal.viscous_strain.first);
    primal.axial_stress = E*elastic_strain;

    primal.overflow_counter++;

    // Terminal condition for reverse sweep.
    adjoint.lambda_history = Dual<LambdaAdj>(LambdaAdj(0.0),LambdaAdj(0.0));

    return make_response(primal);
  }

  MaterialForwardResponse remat_backward_update(const MaterialStepInput& input,
                                                MaterialPrimalState& primal) const {
    const Real dt_abs = std::fabs(input.dt);
    const Real A_real = std::exp(-dt_abs/tau);
    Ratio A_rat(A_real);

    // state in "primal" is at n+1 on entry. Recover state at n.
    // ev_n = (ev_{n+1} - (1-A)*eps_{n+1})/A
    Dual<FixedE> dxx(input.strain_np1*(1-A_rat),0.0);
    primal.viscous_strain = primal.viscous_strain - dxx;
    primal.viscous_strain = primal.viscous_strain/A_rat;

    primal.axial_strain = input.strain_n;

    const Real elastic_strain = primal.axial_strain - Real(primal.viscous_strain.first);
    primal.axial_stress = E*elastic_strain;

    primal.overflow_counter--;

    return make_response(primal);
  }

  LocalAdjointSeed make_objective_seed(const MaterialPrimalState& primal_n) const {
    // Built-in objective contribution used in existing examples/tests:
    // f_n = 0.5*sigma_n^2/E, where sigma_n = E*(eps_n - ev_n).
    LocalAdjointSeed seed;

    const Real elastic_strain_n = primal_n.axial_strain - Real(primal_n.viscous_strain.first);
    const Real sigma_n = E*elastic_strain_n;

    seed.bar_sigma_n = sigma_n/E;
    seed.bar_sigma_np1 = 0.0;
    seed.direct_dE = -0.5*(sigma_n/E)*(sigma_n/E);
    return seed;
  }

  MaterialAdjointResponse adjoint_update(const MaterialStepInput& input,
                                         const MaterialPrimalState& primal_n,
                                         MaterialAdjointState& adjoint,
                                         GradientAccumulator& gradients,
                                         const LocalAdjointSeed& seed) const {
    const Real dt_abs = std::fabs(input.dt);
    const Real A_real = std::exp(-dt_abs/tau);
    Ratio A_rat(A_real);

    const Real ev_n = Real(primal_n.viscous_strain.first);
    const Real eps_n = primal_n.axial_strain;
    const Real eps_np1 = input.strain_np1;
    const Real ev_np1 = A_real*ev_n + (1.0-A_real)*eps_np1;

    // bar(ev_{n+1}) from future plus optional sigma_{n+1} seed.
    Dual<LambdaAdj> bar_ev_np1 = adjoint.lambda_history;
    Dual<LambdaAdj> sigma_np1_term(LambdaAdj(-E*seed.bar_sigma_np1),LambdaAdj(0.0));
    bar_ev_np1 = bar_ev_np1 + sigma_np1_term;

    // dJ/dtau contribution from ev_{n+1} = A*ev_n + (1-A)*eps_{n+1}
    // with dA/dtau = A*dt/tau^2.
    const Real bar_ev_np1_real = Real(bar_ev_np1.first);
    const Real dA_dtau = A_real*(dt_abs/(tau*tau));
    const Real grad_dtau_increment = bar_ev_np1_real*(ev_n - eps_np1)*dA_dtau;
    gradients.df_dtau += grad_dtau_increment;

    // dJ/dE contributions from sigma_n and optional sigma_{n+1}, plus explicit objective term.
    const Real elastic_strain_n = eps_n - ev_n;
    const Real elastic_strain_np1 = eps_np1 - ev_np1;
    const Real grad_dE_increment = seed.bar_sigma_n*elastic_strain_n
                                 + seed.bar_sigma_np1*elastic_strain_np1
                                 + seed.direct_dE;
    gradients.df_dE += grad_dE_increment;

    // Backward recurrence for local history adjoint:
    // lambda_n = A*lambda_{n+1} - E*bar_sigma_n
    if constexpr (std::is_same<LambdaAdj,Real>::value) {
      bar_ev_np1 = bar_ev_np1*A_real;
    } else {
      bar_ev_np1 = bar_ev_np1*A_rat;
    }

    Dual<LambdaAdj> sigma_n_term(LambdaAdj(-E*seed.bar_sigma_n),LambdaAdj(0.0));
    adjoint.lambda_history = bar_ev_np1 + sigma_n_term;

    MaterialAdjointResponse response;
    response.bar_strain_n = E*seed.bar_sigma_n;
    response.bar_strain_np1 = E*seed.bar_sigma_np1 + bar_ev_np1_real*(1.0-A_real);
    response.grad_dtau_increment = grad_dtau_increment;
    response.grad_dE_increment = grad_dE_increment;
    return response;
  }

  // Explicit mode-based update entry point.
  // Forward mode updates primal + resets local adjoint terminal condition.
  // RematBackward mode rematerializes primal state n from state n+1.
  // AdjointBackward mode updates only local adjoint and parameter gradients,
  // using a temporary rematerialized primal state.
  void update(Real lambda, Real& psi, Real* state, Real dt, MaterialUpdateMode mode) {
    check_stretch(lambda);
    const Real strain_n = lambda - 1.0;

    MaterialPrimalState primal = load_primal_state(state);
    MaterialAdjointState adjoint = load_adjoint_state(state);
    GradientAccumulator gradients = load_gradients(state);

    MaterialForwardResponse forward_response;

    if (mode == MaterialUpdateMode::Forward) {
      MaterialStepInput input;
      input.mode = mode;
      input.dt = dt;
      input.strain_n = primal.axial_strain;
      input.strain_np1 = strain_n;
      forward_response = forward_update(input,primal,adjoint);
      store_primal_state(primal,state);
      store_adjoint_state(adjoint,state);
      store_gradients(gradients,state);
      psi = forward_response.energy;
      return;
    }

    if (mode == MaterialUpdateMode::RematBackward) {
      MaterialStepInput input;
      input.mode = mode;
      input.dt = dt;
      input.strain_n = strain_n;
      input.strain_np1 = primal.axial_strain;
      forward_response = remat_backward_update(input,primal);
      store_primal_state(primal,state);
      store_adjoint_state(adjoint,state);
      store_gradients(gradients,state);
      psi = forward_response.energy;
      return;
    }

    // mode == MaterialUpdateMode::AdjointBackward
    // For adjoint-only mode, preserve primal state in storage and perform local
    // rematerialization using a temporary primal copy.
    MaterialPrimalState primal_n = primal;
    MaterialStepInput input;
    input.mode = MaterialUpdateMode::RematBackward;
    input.dt = dt;
    input.strain_n = strain_n;
    input.strain_np1 = primal.axial_strain;

    forward_response = remat_backward_update(input,primal_n);

    if (m_enable_adjoint) {
      LocalAdjointSeed seed = consume_objective_seed(primal_n);
      input.mode = MaterialUpdateMode::AdjointBackward;
      adjoint_update(input,primal_n,adjoint,gradients,seed);
    }

    // Do not overwrite primal state in adjoint-only mode.
    store_adjoint_state(adjoint,state);
    store_gradients(gradients,state);
    psi = forward_response.energy;
  }

  // Backward-compatible update interface.
  // dt >= 0.0: forward material update.
  // dt < 0.0: rematerialize one step backward and apply local adjoint update.
  void update(Real lambda, Real &psi, Real* state, Real dt) {
    if (dt >= 0.0) {
      update(lambda,psi,state,dt,MaterialUpdateMode::Forward);
      return;
    }

    check_stretch(lambda);
    const Real strain_n = lambda - 1.0;

    MaterialPrimalState primal = load_primal_state(state);
    MaterialAdjointState adjoint = load_adjoint_state(state);
    GradientAccumulator gradients = load_gradients(state);

    // On entry in reverse mode, primal is at n+1.
    const Real strain_np1 = primal.axial_strain;

    MaterialStepInput remat_input;
    remat_input.mode = MaterialUpdateMode::RematBackward;
    remat_input.dt = dt;
    remat_input.strain_n = strain_n;
    remat_input.strain_np1 = strain_np1;
    MaterialForwardResponse response = remat_backward_update(remat_input,primal);

    if (m_enable_adjoint) {
      LocalAdjointSeed seed = consume_objective_seed(primal);

      MaterialStepInput adjoint_input;
      adjoint_input.mode = MaterialUpdateMode::AdjointBackward;
      adjoint_input.dt = dt;
      adjoint_input.strain_n = strain_n;
      adjoint_input.strain_np1 = strain_np1;
      adjoint_update(adjoint_input,primal,adjoint,gradients,seed);
    }

    store_primal_state(primal,state);
    store_adjoint_state(adjoint,state);
    store_gradients(gradients,state);
    psi = response.energy;
  }

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state) {
    if ((int(state[kOverflowCounter]) == 0) && (overflow_state.size() >= 2)) {
      state[kOverflowCounter] = Real(mat_overflow_limit);
      // Load in reverse order of storage
      state[kDualLambdaAdjoint] = overflow_state.back(); overflow_state.pop_back();
      state[kDualViscousStrain] = overflow_state.back(); overflow_state.pop_back();
    }
  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    if (int(state[kOverflowCounter]) == mat_overflow_limit) {
      state[kOverflowCounter] = Real(0);
      // Store in forward order to match reverse loading
      overflow_state.push_back(state[kDualViscousStrain]);
      overflow_state.push_back(state[kDualLambdaAdjoint]);
      state[kDualViscousStrain] = Real(0.0);
      state[kDualLambdaAdjoint] = Real(0.0);
    }
  }

  Real get_state_variable(Real* state, std::string state_variable_name) {
    if (state_variable_name == "viscous_strain") {
      FixedE temp;
      load_from_Real(state[kViscousStrain],temp);
      return Real(temp);
    } else if (state_variable_name == "df_dtau") {
      return state[kDfDtau];
    } else if (state_variable_name == "df_dE") {
      return state[kDfDE];
    } else {
      return Real(0.0);
    }
  }

  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    FixedE viscous_temp;
    LambdaAdj lambda_temp;

    field_data[0] = state[kAxialStress];
    field_data[1] = state[kAxialStrain];
    load_from_Real(state[kViscousStrain],viscous_temp); field_data[2] = Real(viscous_temp);
    load_from_Real(state[kDualViscousStrain],viscous_temp); field_data[3] = Real(viscous_temp);
    field_data[4] = state[kOverflowCounter];
    field_data[5] = state[kDfDtau];
    field_data[6] = state[kDfDE];
    load_from_Real(state[kLambdaAdjoint],lambda_temp); field_data[7] = Real(lambda_temp);
    load_from_Real(state[kDualLambdaAdjoint],lambda_temp); field_data[8] = Real(lambda_temp);
  }

  bool is_dead(Real*) { return false; }

}; /* UniaxialViscoelasticity */

#endif /* UNIAXIAL_VISCOELASTICITY_H */
