#include "gtest/gtest.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "Parameters.h"
#include "Fixed.h"
#include "Rational.h"
#include "Truss.h"
#include "UniaxialMaterial.h"
#include "UniaxialViscoelasticity.h"
#include <vector>
#include <iostream>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_V = -6;
const int     EXPONENT_U = -4;
typedef Fixed<RADIX,EXPONENT_V> FixedV;
typedef Fixed<RADIX,EXPONENT_U> FixedU;

namespace {

Real prescribed_right_disp(Real time, Real x, Real) {
  if (time <= 0.0) { return 0.0; }
  return 0.05*(x - 0.0);
}

template <typename ProblemT>
void initialize_visco_truss_problem(ProblemT& problem) {
  const int Nnodes = 2;
  const int Ndofs_per_node = 2;
  const int Nelems = 0;
  const int Nnodes_per_elem = 4;
  double coordinates[4] = { 0.0, 0.0,
                            1.0, 0.0 };
  double velocities[4] = { 0.0, 0.0,
                           0.0, 0.0 };
  bool fixity[4] = { true, true,
                     false, true };
  int connectivity[4] = { 0, 0, 0, 0 };
  int truss_connectivity[2] = { 0, 1 };

  Parameters params;
  params["dt_scale_factor"] = 1.0;
  params["mass_damping_factor"] = 0.0;
  params["density"] = 1.0;
  params["youngs_modulus"] = 1.0;
  params["poissons_ratio"] = 0.25;
  params["truss_density"] = 1.0;
  params["truss_youngs_modulus"] = 1.0;
  params["area"] = 1.0;
  params["relaxation_time"] = 0.35;
  params["mat_overflow_limit"] = 1000000;

  problem.initialize(&coordinates[0],&velocities[0],&fixity[0],Nnodes,Ndofs_per_node,
                     &connectivity[0],Nelems,Nnodes_per_elem,params);
  problem.initialize_truss_elements(&truss_connectivity[0],1,params);
  int moving_node[1] = {1};
  problem.define_displacement_bc(&moving_node[0],1,0,prescribed_right_disp);
  problem.initialize_state();
}

template <typename ProblemT>
double extract_truss_field(ProblemT& problem, const std::string& field_name) {
  const int ntruss = problem.get_num_entities("truss");
  const int nfields = problem.get_num_fields("truss");
  std::vector<double> fields(ntruss*nfields,0.0);
  problem.get_fields("truss",fields.data());
  int field_id = -1;
  for (int i=0; i<nfields; i++) {
    if (field_name == problem.get_field_name("truss",i)) {
      field_id = i;
      break;
    }
  }
  if (field_id < 0) { return 0.0; }
  return fields[field_id];
}

} // end anonymous namespace

TEST(test_System, constructors) {
  // Test creation of System object
  System<Element<Material>,Truss<UniaxialMaterial>,FixedV,FixedU,Rational> problem;
} /* TEST(test_Matrix, constructors) */

TEST(test_System, initialize) {
  System<Element<Material>,Truss<UniaxialMaterial>,FixedV,FixedU,Rational> problem;

  // Define test problem geometry
  const int Nnodes = 9;
  const int Ndofs_per_node = 2;
  const int Nelems = 4;
  const int Nnodes_per_elem = 4;
  double coordinates[18] = { 0.0, 0.0,
                             1.0, 0.0,
                             2.0, 0.0,
		             0.0, 1.0,
                             1.0, 1.0,
                             2.0, 1.0,
		 	     0.0, 2.0,
                             1.0, 2.0,
                             2.0, 2.0 };
  double velocities[18] = { 0.0 };
  bool fixity[18] = { false };
  int connectivity[16] = { 0, 1, 4, 3,
                           1, 2, 5, 4,
                           3, 4, 7, 6,
                           4, 5, 8, 7 };

  // Create parameters object
  Parameters params;

  // Define global parameters
  params["dt_scale_factor"]     =  1.e-5;
  params["body_force_y"]        = -1.e-5;
  params["initial_velocity_y"]  = -1.e-5;
  params["mass_damping_factor"] = 1.0e-9;

  // Material parameters
  params["density"]        = 0.5;
  params["youngs_modulus"] = 215.0e+1; // GPa
  params["poissons_ratio"] = 0.28;

  // Test initialization of the problem object
  problem.initialize(&coordinates[0],&velocities[0],&fixity[0],Nnodes,Ndofs_per_node,
  		     &connectivity[0],Nelems,Nnodes_per_elem,
  		     params);

  // Initialize the problem state at time t=0.0
  std::cout << "Initializing state ... " << std::endl;
  problem.initialize_state();

  // Update the time step
  std::cout << "Updating the time step ... " << std::endl;
  int Nsteps = 10;
  double time = 0.0;
  Real dt = 1.0e-7;
  for (int i = 0; i<Nsteps; i++) {
    time = problem.update_state(dt);
  }
  
} /* TEST(test_Matrix, initialize) */

TEST(test_System, explicit_forward_and_remat_api) {
  System<Element<Material>,Truss<UniaxialMaterial>,FixedV,FixedU,Rational> problem;

  const int Nnodes = 9;
  const int Ndofs_per_node = 2;
  const int Nelems = 4;
  const int Nnodes_per_elem = 4;
  double coordinates[18] = { 0.0, 0.0,
                             1.0, 0.0,
                             2.0, 0.0,
                             0.0, 1.0,
                             1.0, 1.0,
                             2.0, 1.0,
                             0.0, 2.0,
                             1.0, 2.0,
                             2.0, 2.0 };
  double velocities[18] = { 0.0 };
  bool fixity[18] = { false };
  int connectivity[16] = { 0, 1, 4, 3,
                           1, 2, 5, 4,
                           3, 4, 7, 6,
                           4, 5, 8, 7 };

  Parameters params;
  params["dt_scale_factor"]     = 1.e-5;
  params["body_force_y"]        = -1.e-5;
  params["initial_velocity_y"]  = -1.e-5;
  params["mass_damping_factor"] = 1.0e-9;
  params["density"]             = 0.5;
  params["youngs_modulus"]      = 215.0e+1;
  params["poissons_ratio"]      = 0.28;

  problem.initialize(&coordinates[0],&velocities[0],&fixity[0],Nnodes,Ndofs_per_node,
                     &connectivity[0],Nelems,Nnodes_per_elem,
                     params);
  problem.initialize_state();

  const Real dt = 1.0e-7;
  const int Nsteps = 8;
  for (int i=0; i<Nsteps; i++) {
    problem.step_forward(dt);
  }
  EXPECT_GT(problem.get_time(),0.0);

  for (int i=0; i<Nsteps; i++) {
    problem.step_remat_backward();
  }
  EXPECT_NEAR(problem.get_time(),0.0,1.0e-14);

  // No truss elements were defined, so adjoint support should be unsupported.
  EXPECT_EQ(problem.get_truss_adjoint_support(),int(AdjointSupportLevel::Unsupported));
  EXPECT_STREQ(problem.get_truss_adjoint_status(),"unsupported");
}

TEST(test_System, explicit_adjoint_backward_matches_legacy_truss_gradient) {
  using ViscoTrussProblem = System<Element<Material>,
                                   Truss<UniaxialViscoelasticity<Real,Real,Real> >,
                                   Real,Real,Real>;

  ViscoTrussProblem legacy_problem;
  ViscoTrussProblem explicit_problem;
  initialize_visco_truss_problem(legacy_problem);
  initialize_visco_truss_problem(explicit_problem);

  const Real dt = 1.0e-3;
  const int Nsteps = 80;
  for (int i=0; i<Nsteps; i++) {
    legacy_problem.step_forward(dt);
    explicit_problem.step_forward(dt);
  }

  for (int i=0; i<Nsteps; i++) {
    legacy_problem.update_state(-dt);
    explicit_problem.step_adjoint_backward();
  }

  const double legacy_df_dtau = extract_truss_field(legacy_problem,"df_dtau");
  const double explicit_df_dtau = extract_truss_field(explicit_problem,"df_dtau");
  const double legacy_df_dE = extract_truss_field(legacy_problem,"df_dE");
  const double explicit_df_dE = extract_truss_field(explicit_problem,"df_dE");

  const double tol_dtau = 1.0e-10;
  const double tol_dE = 1.0e-10;
  EXPECT_NEAR(legacy_df_dtau,explicit_df_dtau,tol_dtau);
  EXPECT_NEAR(legacy_df_dE,explicit_df_dE,tol_dE);
}
