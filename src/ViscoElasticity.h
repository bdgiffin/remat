#ifndef VISCOELASTICITY_H
#define VISCOELASTICITY_H

#include <math.h>
#include <iostream>
#include <stdlib.h> // exit
#include "Parameters.h"
#include "Dual.h"
#include "Fixed.h"
#include "PassPhase.h"
#include "Rational.h"
#include <limits>
#include <array>
#include <algorithm>
#include <type_traits>
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
  Real adjoint_material_objective_weight = 1.0;
 private:
  enum StateIndex {
    STRESS_XX = 0,
    STRESS_YY = 1,
    STRESS_ZZ = 2,
    STRESS_YZ = 3,
    STRESS_ZX = 4,
    STRESS_XY = 5,
    STRAIN_XX = 6,
    STRAIN_YY = 7,
    STRAIN_XY = 8,
    VISCOUS_XX = 9,
    VISCOUS_YY = 10,
    VISCOUS_XY = 11,
    DUAL_VISCOUS_XX = 12,
    DUAL_VISCOUS_YY = 13,
    DUAL_VISCOUS_XY = 14,
    STIFFNESS_SCALING = 15,
    OVERFLOW_COUNTER = 16,
    DPARAM_RELAXATION_TIME = 17,
    DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT = 18,
    LAMBDA_Q_XX = 19,
    LAMBDA_Q_YY = 20,
    LAMBDA_Q_XY = 21,
    STRESS_SEED_XX = 22,
    STRESS_SEED_YY = 23,
    STRESS_SEED_XY = 24,
    DPARAM_STIFFNESS_SCALING_FACTOR = 25,
    DPARAM_STIFFNESS_SCALING_FACTOR_LOCAL_DEBUG = 26,
    DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG = 27,
    RELAXATION_TIME_LOCAL = 28
  };


  // Deviatoric projector in 2D
  inline void dev2(Real e[3], Real out[3]) {
    Real tr2 = e[0] + e[1];
    out[0] = e[0] - 0.5 * tr2; // eps_xx - tr/2
    out[1] = e[1] - 0.5 * tr2; // eps_yy - tr/2
    out[2] = e[2];             // gxy unchanged
  }

  LongInteger dual_mantissa_abs(const Real&) const { return 0; }

  template<typename T>
  LongInteger dual_mantissa_abs(const T& value) const {
    return std::llabs(LongInteger(value.mantissa));
  }

  void maybe_warn_dual_overflow(FixedE dual_xx, FixedE dual_yy, FixedE dual_xy, PassPhase phase, int overflow_counter) {
    if (!dual_overflow_warn_enable) { return; }
    if (dual_overflow_warn_limit == 0) { return; }
    if (dual_overflow_warn_count >= dual_overflow_warn_limit) { return; }

    const LongInteger max_integer = LongInteger(std::numeric_limits<Integer>::max());
    const LongInteger threshold = LongInteger(dual_overflow_warn_fraction*Real(max_integer));
    const LongInteger abs_xx = dual_mantissa_abs(dual_xx);
    const LongInteger abs_yy = dual_mantissa_abs(dual_yy);
    const LongInteger abs_xy = dual_mantissa_abs(dual_xy);
    const LongInteger max_abs = std::max({abs_xx, abs_yy, abs_xy});

    if (max_abs >= threshold) {
      dual_overflow_warn_count++;
      std::cout
        << "WARNING: ViscoElasticity dual viscous state near overflow"
        << " in phase=" << pass_phase_name(phase)
        << ", overflow_counter=" << overflow_counter
        << ", |mantissa|max=" << max_abs
        << ", threshold=" << threshold
        << ", int32_max=" << max_integer
        << std::endl;
      if (dual_overflow_warn_count == dual_overflow_warn_limit) {
        std::cout << "Further ViscoElasticity dual-overflow warnings suppressed for this run." << std::endl;
      }
    }
  }

  int  mat_overflow_limit = std::numeric_limits<int>::max();
  bool dual_overflow_warn_enable = false;
  Real dual_overflow_warn_fraction = 0.95;
  int dual_overflow_warn_limit = 20;
  int dual_overflow_warn_count = 0;

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

    // 2D relations
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
        std::cout << "Missing viscous parameters: relaxation_time" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (tau<=0.0) {
      std::cout << "ViscoElasticity: need a positive relaxation_time" << std::endl;
      exit(EXIT_FAILURE);
    }

    if (params.count("mat_overflow_limit") > 0) {

      mat_overflow_limit = int(params["mat_overflow_limit"]);
    }
    if (params.count("adjoint_material_objective_weight") > 0) {
      adjoint_material_objective_weight = params["adjoint_material_objective_weight"];
    }
    if (params.count("dual_overflow_warn_enable") > 0) {
      dual_overflow_warn_enable = (params["dual_overflow_warn_enable"] != 0.0);
    }
    if (params.count("dual_overflow_warn_fraction") > 0) {
      dual_overflow_warn_fraction = std::max(0.0,std::min(1.0,params["dual_overflow_warn_fraction"]));
    }
    if (params.count("dual_overflow_warn_limit") > 0) {
      dual_overflow_warn_limit = std::max(0,int(params["dual_overflow_warn_limit"]));
    }
    dual_overflow_warn_count = 0;
  }
    
  // Return the number of state variables for allocation purposes
 int num_state_vars(void) { return 29; }

  // Return the names of all fields
  std::vector<std::string> get_field_names(void) {
    return std::vector<std::string>({
      "stress_xx","stress_yy","stress_zz","stress_yz","stress_zx","stress_xy",
      "strain_xx","strain_yy","strain_xy",
      "viscous_strain_xx","viscous_strain_yy","viscous_strain_xy",
      "viscous_strain_xx_dual","viscous_strain_yy_dual","viscous_strain_xy_dual",
      "stiffness_scaling_factor",
      "overflow_counter",
      "dparam_relaxation_time",
      "dparam_shear_modulus_maxwell_element",
      "lambda_q_xx",
      "lambda_q_yy",
      "lambda_q_xy",
      "adjoint_stress_seed_xx",
      "adjoint_stress_seed_yy",
      "adjoint_stress_seed_xy",
      "dparam_stiffness_scaling_factor",
      "dparam_stiffness_scaling_factor_from_update",
      "dparam_stiffness_scaling_factor_from_direct_seed",
      "relaxation_time_local",
    });
  }

  // Return the density
  Real density(void) { return rho; }
    
  // Initialize the material state
  void initialize(Real* state) {
    state[STRESS_XX] = 0.0;
    state[STRESS_YY] = 0.0;
    state[STRESS_ZZ] = 0.0;
    state[STRESS_YZ] = 0.0;
    state[STRESS_ZX] = 0.0;
    state[STRESS_XY] = 0.0;
    state[STRAIN_XX] = 0.0;
    state[STRAIN_YY] = 0.0;
    state[STRAIN_XY] = 0.0;
    // store primal viscous strains as Fixed_E
    save_as_Real(FixedE(0.0), state[VISCOUS_XX]);
    save_as_Real(FixedE(0.0), state[VISCOUS_YY]);
    save_as_Real(FixedE(0.0), state[VISCOUS_XY]);
    // store dual viscous strains as Fixed_E
    save_as_Real(FixedE(0.0), state[DUAL_VISCOUS_XX]);
    save_as_Real(FixedE(0.0), state[DUAL_VISCOUS_YY]);
    save_as_Real(FixedE(0.0), state[DUAL_VISCOUS_XY]);
    state[STIFFNESS_SCALING] = 1.0;
    state[OVERFLOW_COUNTER] = Real(0);
    state[DPARAM_RELAXATION_TIME] = 0.0;
    state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT] = 0.0;
    state[LAMBDA_Q_XX] = 0.0;
    state[LAMBDA_Q_YY] = 0.0;
    state[LAMBDA_Q_XY] = 0.0;
    state[STRESS_SEED_XX] = 0.0;
    state[STRESS_SEED_YY] = 0.0;
    state[STRESS_SEED_XY] = 0.0;
    state[DPARAM_STIFFNESS_SCALING_FACTOR] = 0.0;
    state[DPARAM_STIFFNESS_SCALING_FACTOR_LOCAL_DEBUG] = 0.0;
    state[DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG] = 0.0;
    state[RELAXATION_TIME_LOCAL] = tau;
  } // initialize()

  // Do we have to keep this?
  // Initialize variable material properties
  void initialize_variable_properties(Real (&x)[2], Real* state, double (*function_xy)(double,double)) {
    // Assign variable stiffness_scaling_factor as a function of initial spatial (x,y) coordinates
    state[STIFFNESS_SCALING] = function_xy(x[0],x[1]);
  } // initialize_variable_properties()

  // Initialize variable relaxation-time values
  void initialize_variable_relaxation_time(Real (&x)[2], Real* state, double (*function_xy)(double,double)) {
    const Real local_tau = function_xy(x[0],x[1]);
    if (local_tau <= 0.0) {
      std::cout << "ViscoElasticity: initialize_variable_relaxation_time requires positive tau, got "
                << local_tau << std::endl;
      exit(EXIT_FAILURE);
    }
    state[RELAXATION_TIME_LOCAL] = local_tau;
  } // initialize_variable_relaxation_time()

  // Backward compatibility path for older materials that used signed dt.
  void update(Real (&F)[2][2],Real &psi, Real* state, Real dt) {
    PassPhase phase = (dt >= 0.0) ? PassPhase::Forward : PassPhase::Backward;
    update(F,psi,state,std::fabs(dt),phase);
  }

  // Update the material state using the current deformation gradient F
  void update(Real (&F)[2][2], Real &psi, Real* state, Real dt, PassPhase phase) {
    if (dt < 0.0) {
      std::cout << "ViscoElasticity::update requires dt >= 0.0 for explicit phases" << std::endl;
      exit(EXIT_FAILURE);
    }

    Real J = F[0][0]*F[1][1] - F[0][1]*F[1][0];
    if (J <= 0.0) {
      std::cout << "F = [[" << F[0][0] << ", " << F[0][1] << "], "<< "[" << F[1][0] << ", " << F[1][1] << "]]\n";
      std::cout << "ERROR: negative Jacobian" << std::endl;
      exit(EXIT_FAILURE);
    }

    const bool reverse_phase = is_reverse_phase(phase);
    const bool run_adjoint = (phase == PassPhase::BackwardAdjoint);
    const Real dt_abs = std::fabs(dt);

    // 2D small strain (engineering shear)
    Real strain_xx = F[0][0] - 1.0;
    Real strain_yy = F[1][1] - 1.0;
    Real strain_xy = F[0][1] + F[1][0];

    // Viscous strain state
    FixedE vs_xx_p, vs_xx_d;
    FixedE vs_yy_p, vs_yy_d;
    FixedE vs_xy_p, vs_xy_d;
    load_from_Real(state[VISCOUS_XX],      vs_xx_p);
    load_from_Real(state[DUAL_VISCOUS_XX], vs_xx_d);
    load_from_Real(state[VISCOUS_YY],      vs_yy_p);
    load_from_Real(state[DUAL_VISCOUS_YY], vs_yy_d);
    load_from_Real(state[VISCOUS_XY],      vs_xy_p);
    load_from_Real(state[DUAL_VISCOUS_XY], vs_xy_d);

    Dual<FixedE> viscous_strain_xx(vs_xx_p, vs_xx_d);
    Dual<FixedE> viscous_strain_yy(vs_yy_p, vs_yy_d);
    Dual<FixedE> viscous_strain_xy(vs_xy_p, vs_xy_d);

    int overflow_counter = int(state[OVERFLOW_COUNTER]);
    maybe_warn_dual_overflow(vs_xx_d,vs_yy_d,vs_xy_d,phase,overflow_counter);
    Real previous_strain[3] = { state[STRAIN_XX], state[STRAIN_YY], state[STRAIN_XY] };
    Real dev_previous_strain[3] = { 0.0, 0.0, 0.0 };
    Real dev_strain_np1[3] = { 0.0, 0.0, 0.0 };
    Real lambda_q_np1[3] = { state[LAMBDA_Q_XX], state[LAMBDA_Q_YY], state[LAMBDA_Q_XY] };
    Real dparam_tau = state[DPARAM_RELAXATION_TIME];
    Real dparam_mu_e = state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT];
    Real dparam_stiffness_scaling = state[DPARAM_STIFFNESS_SCALING_FACTOR];
    Real dparam_stiffness_scaling_local_debug = state[DPARAM_STIFFNESS_SCALING_FACTOR_LOCAL_DEBUG];
    Real dparam_stiffness_scaling_direct_debug = state[DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG];
    const Real stress_seed_xx = state[STRESS_SEED_XX];
    const Real stress_seed_yy = state[STRESS_SEED_YY];
    const Real stress_seed_xy = state[STRESS_SEED_XY];

    const Real tau_local = state[RELAXATION_TIME_LOCAL];
    if (tau_local <= 0.0) {
      std::cout << "ViscoElasticity: local relaxation_time must remain positive, got "
                << tau_local << std::endl;
      exit(EXIT_FAILURE);
    }
    const Real A = std::exp(-dt_abs/tau_local);
    Ratio A_rat(A);

    if (phase == PassPhase::Forward) {
      viscous_strain_xx = viscous_strain_xx * A_rat;
      viscous_strain_yy = viscous_strain_yy * A_rat;
      viscous_strain_xy = viscous_strain_xy * A_rat;

      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx(dev_previous_strain[0]*(1-A_rat),0.0);
      Dual<FixedE> dyy(dev_previous_strain[1]*(1-A_rat),0.0);
      Dual<FixedE> dxy(dev_previous_strain[2]*(1-A_rat),0.0);
      viscous_strain_xx = viscous_strain_xx + dxx;
      viscous_strain_yy = viscous_strain_yy + dyy;
      viscous_strain_xy = viscous_strain_xy + dxy;

      // Forward sweep terminal state for reverse adjoint.
      state[LAMBDA_Q_XX] = 0.0;
      state[LAMBDA_Q_YY] = 0.0;
      state[LAMBDA_Q_XY] = 0.0;

      overflow_counter++;
    } else if (reverse_phase) {
      if (run_adjoint) {
        dev2(previous_strain,dev_strain_np1); // epsilon_{n+1} deviatoric part
      }

      // Inverse rematerialization of viscous history.
      dev2(previous_strain,dev_previous_strain);
      Dual<FixedE> dxx(dev_previous_strain[0]*(1-A_rat),0.0);
      Dual<FixedE> dyy(dev_previous_strain[1]*(1-A_rat),0.0);
      Dual<FixedE> dxy(dev_previous_strain[2]*(1-A_rat),0.0);
      viscous_strain_xx = viscous_strain_xx - dxx;
      viscous_strain_yy = viscous_strain_yy - dyy;
      viscous_strain_xy = viscous_strain_xy - dxy;

      previous_strain[0] = strain_xx;
      previous_strain[1] = strain_yy;
      previous_strain[2] = strain_xy;

      viscous_strain_xx = viscous_strain_xx / A_rat;
      viscous_strain_yy = viscous_strain_yy / A_rat;
      viscous_strain_xy = viscous_strain_xy / A_rat;

      if (run_adjoint) {
        const Real q_n[3] = { Real(viscous_strain_xx.first), Real(viscous_strain_yy.first), Real(viscous_strain_xy.first) };
        Real dev_strain_n[3] = { 0.0, 0.0, 0.0 };
        dev2(previous_strain,dev_strain_n);
        const Real dev_elastic_n[3] = {
          dev_strain_n[0] - q_n[0],
          dev_strain_n[1] - q_n[1],
          dev_strain_n[2] - q_n[2]
        };

        const Real stiffness_scaling_factor = state[STIFFNESS_SCALING];
        const Real mu2_scaled = stiffness_scaling_factor * mu2;
        const Real mu_scaled  = stiffness_scaling_factor * mu;
        const Real lam_scaled = stiffness_scaling_factor * lam;
        const Real mu2_e_scaled = stiffness_scaling_factor * mu2_e;
        const Real mu_e_scaled  = stiffness_scaling_factor * mu_e;

        const Real stress_xx_eq_n = (lam_scaled + mu2_scaled)*previous_strain[0] + lam_scaled*previous_strain[1];
        const Real stress_yy_eq_n = lam_scaled*previous_strain[0] + (lam_scaled + mu2_scaled)*previous_strain[1];
        const Real stress_xy_eq_n = mu_scaled*previous_strain[2];
        const Real stress_xx_n = stress_xx_eq_n + mu2_e_scaled*dev_elastic_n[0];
        const Real stress_yy_n = stress_yy_eq_n + mu2_e_scaled*dev_elastic_n[1];
        const Real stress_xy_n = stress_xy_eq_n + mu_e_scaled*dev_elastic_n[2];

        // Built-in objective seed: 0.5*(sxx^2 + syy^2 + 2*sxy^2), plus external/global stress seeds.
        const Real bar_sigma_xx = adjoint_material_objective_weight*stress_xx_n + stress_seed_xx;
        const Real bar_sigma_yy = adjoint_material_objective_weight*stress_yy_n + stress_seed_yy;
        const Real bar_sigma_xy = 2.0*adjoint_material_objective_weight*stress_xy_n + stress_seed_xy;

        dparam_tau += lambda_q_np1[0]*(q_n[0] - dev_strain_np1[0])*A*(dt_abs/(tau_local*tau_local));
        dparam_tau += lambda_q_np1[1]*(q_n[1] - dev_strain_np1[1])*A*(dt_abs/(tau_local*tau_local));
        dparam_tau += lambda_q_np1[2]*(q_n[2] - dev_strain_np1[2])*A*(dt_abs/(tau_local*tau_local));

        const Real dsigma_xx_dscale =
          (lam + mu2)*previous_strain[0] + lam*previous_strain[1] + mu2_e*dev_elastic_n[0];
        const Real dsigma_yy_dscale =
          lam*previous_strain[0] + (lam + mu2)*previous_strain[1] + mu2_e*dev_elastic_n[1];
        const Real dsigma_xy_dscale =
          mu*previous_strain[2] + mu_e*dev_elastic_n[2];
        const Real dparam_stiffness_scaling_local_increment =
          bar_sigma_xx*dsigma_xx_dscale +
          bar_sigma_yy*dsigma_yy_dscale +
          bar_sigma_xy*dsigma_xy_dscale;
        dparam_stiffness_scaling += dparam_stiffness_scaling_local_increment;
        dparam_stiffness_scaling_local_debug += dparam_stiffness_scaling_local_increment;

        dparam_mu_e += bar_sigma_xx*(2.0*stiffness_scaling_factor*dev_elastic_n[0]);
        dparam_mu_e += bar_sigma_yy*(2.0*stiffness_scaling_factor*dev_elastic_n[1]);
        dparam_mu_e += bar_sigma_xy*(stiffness_scaling_factor*dev_elastic_n[2]);

        state[LAMBDA_Q_XX] = A*lambda_q_np1[0] - mu2_e_scaled*bar_sigma_xx;
        state[LAMBDA_Q_YY] = A*lambda_q_np1[1] - mu2_e_scaled*bar_sigma_yy;
        state[LAMBDA_Q_XY] = A*lambda_q_np1[2] - mu_e_scaled*bar_sigma_xy;
      }

      overflow_counter--;
    } else {
      std::cout << "ViscoElasticity::update received unsupported phase: "
                << pass_phase_name(phase) << std::endl;
      exit(EXIT_FAILURE);
    }

    // Stress update from rematerialized state.
    const Real stiffness_scaling_factor = state[STIFFNESS_SCALING];
    const Real mu2_scaled = stiffness_scaling_factor * mu2;
    const Real mu_scaled  = stiffness_scaling_factor * mu;
    const Real lam_scaled = stiffness_scaling_factor * lam;
    const Real mu2_e_scaled = stiffness_scaling_factor * mu2_e;
    const Real mu_e_scaled  = stiffness_scaling_factor * mu_e;

    const Real strain_now_xx = previous_strain[0];
    const Real strain_now_yy = previous_strain[1];
    const Real strain_now_xy = previous_strain[2];

    const Real stress_xx_eq = (lam_scaled + mu2_scaled)*strain_now_xx + lam_scaled*strain_now_yy;
    const Real stress_yy_eq = lam_scaled*strain_now_xx + (lam_scaled + mu2_scaled)*strain_now_yy;
    const Real stress_xy_eq = mu_scaled * strain_now_xy;

    dev2(previous_strain,dev_previous_strain);
    const Real dev_elastic_strain_xx = dev_previous_strain[0] - Real(viscous_strain_xx.first);
    const Real dev_elastic_strain_yy = dev_previous_strain[1] - Real(viscous_strain_yy.first);
    const Real dev_elastic_strain_xy = dev_previous_strain[2] - Real(viscous_strain_xy.first);

    const Real stress_xx_Maxwell = mu2_e_scaled*dev_elastic_strain_xx;
    const Real stress_yy_Maxwell = mu2_e_scaled*dev_elastic_strain_yy;
    const Real stress_xy_Maxwell = mu_e_scaled*dev_elastic_strain_xy;

    const Real stress_xx = stress_xx_eq + stress_xx_Maxwell;
    const Real stress_yy = stress_yy_eq + stress_yy_Maxwell;
    const Real stress_xy = stress_xy_eq + stress_xy_Maxwell;

    // Check the updated dual viscous history as well, so we can warn if a
    // value grows dangerously between overflow checkpoint stores.
    maybe_warn_dual_overflow(
      viscous_strain_xx.second,
      viscous_strain_yy.second,
      viscous_strain_xy.second,
      phase,
      overflow_counter
    );

    state[STRESS_ZZ] = 0.0;
    state[STRESS_YZ] = 0.0;
    state[STRESS_ZX] = 0.0;
    state[STRESS_XX] = stress_xx;
    state[STRESS_YY] = stress_yy;
    state[STRESS_XY] = stress_xy;
    state[STRAIN_XX] = previous_strain[0];
    state[STRAIN_YY] = previous_strain[1];
    state[STRAIN_XY] = previous_strain[2];

    save_as_Real(viscous_strain_xx.first,  state[VISCOUS_XX]);
    save_as_Real(viscous_strain_yy.first,  state[VISCOUS_YY]);
    save_as_Real(viscous_strain_xy.first,  state[VISCOUS_XY]);
    save_as_Real(viscous_strain_xx.second, state[DUAL_VISCOUS_XX]);
    save_as_Real(viscous_strain_yy.second, state[DUAL_VISCOUS_YY]);
    save_as_Real(viscous_strain_xy.second, state[DUAL_VISCOUS_XY]);

    state[OVERFLOW_COUNTER] = Real(overflow_counter);
    state[DPARAM_RELAXATION_TIME] = dparam_tau;
    state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT] = dparam_mu_e;
    state[DPARAM_STIFFNESS_SCALING_FACTOR] = dparam_stiffness_scaling;
    state[DPARAM_STIFFNESS_SCALING_FACTOR_LOCAL_DEBUG] = dparam_stiffness_scaling_local_debug;
    state[DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG] = dparam_stiffness_scaling_direct_debug;
    state[STRESS_SEED_XX] = 0.0;
    state[STRESS_SEED_YY] = 0.0;
    state[STRESS_SEED_XY] = 0.0;

    const Real psi_eq = 0.5*(stress_xx_eq*strain_now_xx + stress_yy_eq*strain_now_yy + stress_xy_eq*strain_now_xy);
    const Real psi_Maxwell = 0.5*(stress_xx_Maxwell*dev_elastic_strain_xx +
                                  stress_yy_Maxwell*dev_elastic_strain_yy +
                                  stress_xy_Maxwell*dev_elastic_strain_xy);
    psi = psi_eq + psi_Maxwell;
  } // update()

  // Conditionally load material history parameters from memory
  void load_state(Real* state, std::vector<Real>& overflow_state)  {
    if ((int(state[OVERFLOW_COUNTER]) == 0) && (overflow_state.size() > 0)) {
      state[OVERFLOW_COUNTER] = Real(mat_overflow_limit);
      // load in reverse order
      state[DUAL_VISCOUS_XY] = overflow_state.back(); overflow_state.pop_back();
      state[DUAL_VISCOUS_YY] = overflow_state.back(); overflow_state.pop_back();
      state[DUAL_VISCOUS_XX] = overflow_state.back(); overflow_state.pop_back();
    }
  }

  // Conditionally store material history parameters in memory
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    if (int(state[OVERFLOW_COUNTER]) == mat_overflow_limit) {
      state[OVERFLOW_COUNTER] = Real(0);
      // store in forward order
      overflow_state.push_back(state[DUAL_VISCOUS_XX]);
      overflow_state.push_back(state[DUAL_VISCOUS_YY]);
      overflow_state.push_back(state[DUAL_VISCOUS_XY]);
      state[DUAL_VISCOUS_XX] = Real(0.0);
      state[DUAL_VISCOUS_YY] = Real(0.0);
      state[DUAL_VISCOUS_XY] = Real(0.0);
    }
  }

  int adjoint_num_params(void) const { return 3; }

  const char* adjoint_param_name(int i) const {
    switch (i) {
      case 0: return "relaxation_time";
      case 1: return "shear_modulus_Maxwell_element";
      case 2: return "stiffness_scaling_factor";
      default: return "";
    }
  }

  void adjoint_add_stress_seed(Real* state, Real seed_xx, Real seed_yy, Real seed_xy) {
    state[STRESS_SEED_XX] += seed_xx;
    state[STRESS_SEED_YY] += seed_yy;
    state[STRESS_SEED_XY] += seed_xy;
  }

  // Second-kick stress seeds act at the current rematerialized state (n+1),
  // so their history contribution must enter lambda_q_(n+1) directly.
  void adjoint_add_history_seed_from_stress(Real* state,
                                            Real seed_xx, Real seed_yy, Real seed_xy) {
    const Real stiffness_scaling_factor = state[STIFFNESS_SCALING];
    const Real mu2_e_scaled = stiffness_scaling_factor * mu2_e;
    const Real mu_e_scaled  = stiffness_scaling_factor * mu_e;
    state[LAMBDA_Q_XX] += -mu2_e_scaled*seed_xx;
    state[LAMBDA_Q_YY] += -mu2_e_scaled*seed_yy;
    state[LAMBDA_Q_XY] += -mu_e_scaled*seed_xy;
  }

  void adjoint_objective_seed(Real*) { }

  Real adjoint_get_param_gradient(const Real* state, int i) const {
    if (i == 0) { return state[DPARAM_RELAXATION_TIME]; }
    if (i == 1) { return state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT]; }
    if (i == 2) { return state[DPARAM_STIFFNESS_SCALING_FACTOR]; }
    return Real(0.0);
  }

  void adjoint_clear_step_seed(Real* state) {
    state[STRESS_SEED_XX] = 0.0;
    state[STRESS_SEED_YY] = 0.0;
    state[STRESS_SEED_XY] = 0.0;
  }

  void adjoint_pullback_stress_to_strain(const Real* state,
                                         Real stress_seed_xx, Real stress_seed_yy, Real stress_seed_xy,
                                         Real& strain_seed_xx, Real& strain_seed_yy, Real& strain_seed_xy) const {
    const Real stiffness_scaling_factor = state[STIFFNESS_SCALING];
    const Real mu2_scaled = stiffness_scaling_factor * mu2;
    const Real mu_scaled  = stiffness_scaling_factor * mu;
    const Real lam_scaled = stiffness_scaling_factor * lam;
    const Real mu_e_scaled = stiffness_scaling_factor * mu_e;

    const Real dxx_dexx = lam_scaled + mu2_scaled + mu_e_scaled;
    const Real dxx_deyy = lam_scaled - mu_e_scaled;
    const Real dyy_dexx = lam_scaled - mu_e_scaled;
    const Real dyy_deyy = lam_scaled + mu2_scaled + mu_e_scaled;
    const Real dxy_dgxy = mu_scaled + mu_e_scaled;

    strain_seed_xx += dxx_dexx*stress_seed_xx + dyy_dexx*stress_seed_yy;
    strain_seed_yy += dxx_deyy*stress_seed_xx + dyy_deyy*stress_seed_yy;
    strain_seed_xy += dxy_dgxy*stress_seed_xy;
  }

  void adjoint_add_direct_param_seed_from_stress(Real* state,
                                                  Real stress_seed_xx, Real stress_seed_yy, Real stress_seed_xy) {
    FixedE qxx_fix, qyy_fix, qxy_fix;
    load_from_Real(state[VISCOUS_XX], qxx_fix);
    load_from_Real(state[VISCOUS_YY], qyy_fix);
    load_from_Real(state[VISCOUS_XY], qxy_fix);

    const Real exx = state[STRAIN_XX];
    const Real eyy = state[STRAIN_YY];
    const Real gxy = state[STRAIN_XY];
    const Real tr2 = exx + eyy;
    const Real dev_xx = exx - 0.5*tr2;
    const Real dev_yy = eyy - 0.5*tr2;
    const Real dev_xy = gxy;

    const Real dev_elastic_xx = dev_xx - Real(qxx_fix);
    const Real dev_elastic_yy = dev_yy - Real(qyy_fix);
    const Real dev_elastic_xy = dev_xy - Real(qxy_fix);

    const Real stiffness_scaling_factor = state[STIFFNESS_SCALING];
    state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT] +=
      stress_seed_xx*(2.0*stiffness_scaling_factor*dev_elastic_xx) +
      stress_seed_yy*(2.0*stiffness_scaling_factor*dev_elastic_yy) +
      stress_seed_xy*(stiffness_scaling_factor*dev_elastic_xy);

    const Real dsigma_xx_dscale = (lam + mu2)*exx + lam*eyy + mu2_e*dev_elastic_xx;
    const Real dsigma_yy_dscale = lam*exx + (lam + mu2)*eyy + mu2_e*dev_elastic_yy;
    const Real dsigma_xy_dscale = mu*gxy + mu_e*dev_elastic_xy;
    const Real dparam_stiffness_scaling_direct_increment =
      stress_seed_xx*dsigma_xx_dscale +
      stress_seed_yy*dsigma_yy_dscale +
      stress_seed_xy*dsigma_xy_dscale;
    state[DPARAM_STIFFNESS_SCALING_FACTOR] += dparam_stiffness_scaling_direct_increment;
    state[DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG] += dparam_stiffness_scaling_direct_increment;
  }


  // Copy state variable data to field data
  void get_fields(Real* state, double* field_data) {
    field_data[0]  = state[STRESS_XX];
    field_data[1]  = state[STRESS_YY];
    field_data[2]  = state[STRESS_ZZ];
    field_data[3]  = state[STRESS_YZ];
    field_data[4]  = state[STRESS_ZX];
    field_data[5]  = state[STRESS_XY];
    field_data[6]  = state[STRAIN_XX];
    field_data[7]  = state[STRAIN_YY];
    field_data[8]  = state[STRAIN_XY];
    FixedE temp;
    load_from_Real(state[VISCOUS_XX]     ,temp); field_data[9]  = Real(temp);
    load_from_Real(state[VISCOUS_YY]     ,temp); field_data[10] = Real(temp);
    load_from_Real(state[VISCOUS_XY]     ,temp); field_data[11] = Real(temp);
    load_from_Real(state[DUAL_VISCOUS_XX],temp); field_data[12] = Real(temp);
    load_from_Real(state[DUAL_VISCOUS_YY],temp); field_data[13] = Real(temp);
    load_from_Real(state[DUAL_VISCOUS_XY],temp); field_data[14] = Real(temp);
    field_data[15] = state[STIFFNESS_SCALING];
    field_data[16] = state[OVERFLOW_COUNTER];
    field_data[17] = state[DPARAM_RELAXATION_TIME];
    field_data[18] = state[DPARAM_SHEAR_MODULUS_MAXWELL_ELEMENT];
    field_data[19] = state[LAMBDA_Q_XX];
    field_data[20] = state[LAMBDA_Q_YY];
    field_data[21] = state[LAMBDA_Q_XY];
    field_data[22] = state[STRESS_SEED_XX];
    field_data[23] = state[STRESS_SEED_YY];
    field_data[24] = state[STRESS_SEED_XY];
    field_data[25] = state[DPARAM_STIFFNESS_SCALING_FACTOR];
    field_data[26] = state[DPARAM_STIFFNESS_SCALING_FACTOR_LOCAL_DEBUG];
    field_data[27] = state[DPARAM_STIFFNESS_SCALING_FACTOR_DIRECT_DEBUG];
    field_data[28] = state[RELAXATION_TIME_LOCAL];
  }

  // Return the initial sound speed
  //Real initial_sound_speed(void) { return sqrt(pmod/rho); }

}; /* ViscoElasticity */

#endif // VISCOELASTICITY_H
