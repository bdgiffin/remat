// ======================================================================== //
// C API functions for interacting with the REMAT library                   //
// ======================================================================== //

#include "types.h"
#include "Parameters.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "Fixed.h"
#include "Rational.h"
#include "System.h"
#include "Truss.h"
#include "UniaxialMaterial.h"
#include <iostream>
#include <stdio.h>

// Declare standard Fixed-precision numbers
const int          RADIX = 10;
const int     EXPONENT_V = -6;
const int     EXPONENT_U = -6;
typedef Fixed<RADIX,EXPONENT_V> FixedV;
typedef Fixed<RADIX,EXPONENT_U> FixedU;

// global parameter list
Parameters params;

// global instance of the system object
System<Element<Material>,Truss<UniaxialMaterial>,Real,Real,Real>         remat_float;
System<Element<Material>,Truss<UniaxialMaterial>,FixedV,FixedU,Rational> remat_fixed;
SystemBase* remat = &remat_float;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {
  
  // ------------------------------------------------------------------------ //

  // Set which integrator type to use
  void set_integrator_type(const char* integrator_type) {
    std::string integrator_type_string(integrator_type);
    if (integrator_type_string == "float") {
      remat = &remat_float;
    } else if (integrator_type_string == "fixed") {
      remat = &remat_fixed;
    } else {
      std::cerr << "Invalid integrator type specified!" << std::endl;
      exit(EXIT_FAILURE);
    }
    
  } // set_integrator_type()
  
  // ------------------------------------------------------------------------ //

  // Define a named parameter and its corresponding value
  //   name: the string-valued name of the prameter being defined
  //  value: the real-valued parameter constant being defined
  void define_parameter(const char* name, double value) {
    params[name] = value;
  } // define_parameter()
  
  // ------------------------------------------------------------------------ //

  // Define the problem geometry and initialize the System object
  void define_geometry(double *coordinates, double *velocities, bool *fixity,
		       int *connectivity, size_t Nnodes, size_t Nelems) {
    int Ndofs_per_node  = 2;
    int Nnodes_per_elem = 4;
    remat->initialize(coordinates, velocities, fixity, Nnodes, Ndofs_per_node,
	 	      connectivity, Nelems, Nnodes_per_elem,
	 	      params);
  } // define_geometry()
  
  // ------------------------------------------------------------------------ //

  // Define a new contact interaction
  void define_truss_elements(int *truss_connectivity, size_t Ntruss) {
    remat->initialize_truss_elements(truss_connectivity,Ntruss,params);
  } // define_truss_elements()
  
  // ------------------------------------------------------------------------ //

  // Define a new contact interaction
  void define_contact_interaction(int* node_ids, int* segment_connectivity,
		                  int Nnodes, int Nsegments) {
    remat->initialize_contact(node_ids,segment_connectivity,Nnodes,Nsegments,params);
  } // define_contact_interaction()
  
  // ------------------------------------------------------------------------ //

  // Define new point masses
  void define_point_mass(int *point_ids, double *point_mass, size_t Npoints) {
    remat->initialize_point_mass(point_ids,point_mass,Npoints,params);
  } // define_point_mass()
  
  // ------------------------------------------------------------------------ //

  // Initialize the System state at the inidicated initial analysis time
  void initialize() {
    remat->initialize_state();
  } // initialize()
  
  // ------------------------------------------------------------------------ //

  // Update the System state
  double update_state(double dt, int Nsub_steps) {
    double time;
    for (int i=0; i<Nsub_steps; i++) {
      time = remat->update_state(dt);
    }
    return time;
  } // update_state()
  
  // ------------------------------------------------------------------------ //

  // Request the System field data
  double get_field_data(double *ux, double *uy,
	 	        double *vx, double *vy,
		        double *fx, double *fy,
			double *dual_ux, double *dual_uy,
	 	        double *dual_vx, double *dual_vy,
		        double *sxx, double *syy, double *sxy, double *pressure) {
    return remat->get_field_data(ux,uy,vx,vy,fx,fy,dual_ux,dual_uy,dual_vx,dual_vy,sxx,syy,sxy,pressure);
  } // get_field_data()
  
  // ------------------------------------------------------------------------ //
  
} // extern "C"

// ======================================================================== //
