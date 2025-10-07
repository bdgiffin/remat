// ======================================================================== //
// C API functions for interacting with the REMAT library                   //
// ======================================================================== //

#include "types.h"
#include "Parameters.h"
#include "System.h"
#include "Element.h"
#include "Material.h"
#include "ViscoElasticity.h"
#include "Fixed.h"
#include "Rational.h"
#include "System.h"
#include "Truss.h"
#include "UniaxialMaterial.h"
#include "UniaxialViscoplasticity.h"
#include <chrono>
#include <iostream>
#include <stdio.h>

// Declare element and material types
typedef Element<Material>                                 ElementT;
typedef Element<ViscoElasticity<Real,Real> >              ElementT_float_visco;
typedef Element<ViscoElasticity<Fixed_E,Rational> >       ElementT_fixed_visco;
typedef Truss<UniaxialViscoplasticity<Real,Real> >        TrussT_float;
typedef Truss<UniaxialViscoplasticity<Fixed_E,Rational> > TrussT_fixed;

// global parameter list
Parameters params;

// global instance of the system object
System<ElementT,TrussT_float,Real,Real,Real>           remat_float;
System<ElementT,TrussT_fixed,Fixed_V,Fixed_U,Rational> remat_fixed;
System<ElementT,TrussT_fixed,Fixed_V,Fixed_U,Real>     remat_mixed;
System<ElementT_float_visco,TrussT_float,Real,Real,Real>           remat_float_visco;
System<ElementT_fixed_visco,TrussT_fixed,Fixed_V,Fixed_U,Rational> remat_fixed_visco;
SystemBase* remat = &remat_float;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {
  
  // ------------------------------------------------------------------------ //

  // Set which integrator type to use
  void set_integrator_type(const char* integrator_type) {
    std::string integrator_type_string(integrator_type);
    if        (integrator_type_string == "float") {
      remat = &remat_float;
    } else if (integrator_type_string == "fixed") {
      remat = &remat_fixed;
    } else if (integrator_type_string == "mixed") {
      remat = &remat_mixed;
    } else if (integrator_type_string == "float_visco") {
      remat = &remat_mixed;
    } else if (integrator_type_string == "fixed_visco") {
      remat = &remat_mixed;
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

  // Initialize variable material properties
  void initialize_variable_properties(double (*function_xy)(double,double)) {
    remat->initialize_variable_properties(function_xy);
  } // initialize_variable_properties()
  
  // ------------------------------------------------------------------------ //

  // Initialize the System state at the inidicated initial analysis time
  void initialize() {
    remat->initialize_state();
  } // initialize()
  
  // ------------------------------------------------------------------------ //

  // Update the System state
  double update_state(double dt, int Nsub_steps) {
    // Record the starting wall time
    auto start_time = std::chrono::high_resolution_clock::now();
    
    double time;
    for (int i=0; i<Nsub_steps; i++) {
      time = remat->update_state(dt);
    }
    
    // Record the ending wall time
    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed wall time duration
    auto elapsed_duration = end_time - start_time;

    // Convert to milliseconds and get the count
    auto ms_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_duration).count();

    std::cout << "Elapsed time: " << ms_elapsed << " milliseconds" << std::endl;
    
    return time;
    
  } // update_state()
  
  // ------------------------------------------------------------------------ //

  // Request the System field data
  double get_field_data(double *ux, double *uy,
	 	        double *vx, double *vy,
		        double *fx, double *fy,
			double *dual_ux, double *dual_uy,
	 	        double *dual_vx, double *dual_vy,
		        double *sxx, double *syy, double *sxy,
			double *pressure, double *stiffness_scaling_factor,
			double *system_state, double *eqps, bool *is_dead) {
    return remat->get_field_data(ux,uy,vx,vy,fx,fy,dual_ux,dual_uy,dual_vx,dual_vy,
				 sxx,syy,sxy,pressure,stiffness_scaling_factor,system_state,
				 eqps,is_dead);
  } // get_field_data()
  
  // ------------------------------------------------------------------------ //

  // Request the number of spatial dimensions
  int get_num_dim(void) { return remat->get_num_dim(); }
  
  // ------------------------------------------------------------------------ //

  // Request the number of entities of the specified type
  int get_num_entities(const char* entity_type) {
    return remat->get_num_entities(entity_type);
  } // get_num_entities()
  
  // ------------------------------------------------------------------------ //

  // Request the number of fields defined for entities of the specified type
  int get_num_fields(const char* entity_type) {
    return remat->get_num_fields(entity_type);
  } // get_num_fields()
  
  // ------------------------------------------------------------------------ //

  // Request the name of the indicated field ID for entities of the specified type
  const char* get_field_name(const char* entity_type, int field_id) {
    return remat->get_field_name(entity_type,field_id);
  } // get_field_name()
  
  // ------------------------------------------------------------------------ //

  // Request the coordinates of all nodes
  void get_node_coords(double* coords, bool deformed) {
    remat->get_node_coords(coords,deformed);
  } // get_node_coords()
  
  // ------------------------------------------------------------------------ //

  // Request the connectivity data for all entities of the specified type
  void get_connectivity(const char* entity_type, int* connectivity) {
    remat->get_connectivity(entity_type,connectivity);
  } // get_connectivity()
  
  // ------------------------------------------------------------------------ //

  // Request data defining all fields for all entities of the specified type
  void get_fields(const char* entity_type, double* field_data) {
    remat->get_fields(entity_type,field_data);
  } // get_fields()
  
  // ------------------------------------------------------------------------ //

  // Request the current analysis time
  double get_time(void) { return remat->get_time(); }
  
  // ------------------------------------------------------------------------ //
  
} // extern "C"

// ======================================================================== //
