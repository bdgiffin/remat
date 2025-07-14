#include "gtest/gtest.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "Parameters.h"
#include "Fixed.h"
#include "Rational.h"
#include "Truss.h"
#include "UniaxialMaterial.h"
#include <vector>
#include <iostream>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_V = -6;
const int     EXPONENT_U = -4;
typedef Fixed<RADIX,EXPONENT_V> FixedV;
typedef Fixed<RADIX,EXPONENT_U> FixedU;

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
