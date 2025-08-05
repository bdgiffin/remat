#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.h"
#include "ContactInteraction.h"
#include "Dual.h"
#include "Integrator.h"
#include "Parameters.h"
#include <limits>
#include <vector>
#include <algorithm> // For std::fill
#include <iostream>
#include <math.h>

// Declare a non-templated base class to enable instantiation with variably typed 
struct SystemBase {
  
  // default constructor method for a system object
  SystemBase() { }
    
  // ===================================================================== //

  // Procedure to initialize the system
  virtual void initialize(double* new_coordinates, double* new_velocities, bool* new_fixity, int new_Nnodes, int new_Ndofs_per_node,
		          int* new_connectivity, int new_Nelems, int new_Nnodes_per_elem,
		          Parameters& params) = 0;
  
  // ===================================================================== //

  // Procedure to initialize a new contact interaction
  virtual void initialize_contact(int* node_ids, int* segment_connectivity, int new_Nnodes, int new_Nsegments, Parameters& params) = 0;
  
  // ===================================================================== //

  // Procedure to initialize truss elements
  virtual void initialize_truss_elements(int* new_truss_connectivity, int new_Ntruss, Parameters& params) = 0;
  
  // ===================================================================== //

  // Procedure to initialize new point masses
  virtual void initialize_point_mass(int* new_point_ids, double* new_point_mass, int new_Npoints, Parameters& params) = 0;
  
  // ===================================================================== //

  // Procedure to initialize variable material stiffness properties
  virtual void initialize_variable_properties(double (*function_xy)(double,double)) = 0;
  
  // ===================================================================== //

  // Procedure to initialize the system state at the indicated time
  virtual void initialize_state(void) = 0;
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  virtual double update_state(Real dt) = 0;
  
  // ===================================================================== //

  // Procedure to get the current system state data
  virtual double get_field_data(double *ux, double *uy,
	 	                double *vx, double *vy,
		                double *fx, double *fy,
			        double *dual_ux, double *dual_uy,
	 	                double *dual_vx, double *dual_vy,
			        double *sxx, double *syy, double *sxy,
				double *pressure, double *stiffness_scaling_factor,
				double *system_state, double *eqps, bool *is_dead) = 0;
  
  // ===================================================================== //
  
  // Request the number of spatial dimensions
  virtual int get_num_dim(void) = 0;
  
  // ===================================================================== //

  // Request the number of entities of the specified type
  virtual int get_num_entities(std::string entity_type) = 0;
  
  // ===================================================================== //

  // Request the number of fields defined for entities of the specified type
  virtual int get_num_fields(std::string entity_type) = 0;
  
  // ===================================================================== //

  // Request the name of the indicated field ID for entities of the specified type
  virtual const char* get_field_name(std::string entity_type, int field_id) = 0;
  
  // ===================================================================== //
  
  // Request the coordinates of all nodes
  virtual void get_node_coords(double* coords, bool deformed) = 0;
  
  // ===================================================================== //
  
  // Request the connectivity data for all entities of the specified type
  virtual void get_connectivity(std::string entity_type, int* connectivity) = 0;
  
  // ===================================================================== //

  // Request data defining all fields for all entities of the specified type
  virtual void get_fields(std::string entity_type, double* field_data) = 0;
  
  // ===================================================================== //

  // Request the current analysis time
  virtual double get_time(void) = 0;
  
  // ===================================================================== //

}; // SystemBase

// ....................................................................... //

template<class Element_T, class Truss_T, class FixedV, class FixedU, class Ratio>
struct System : public SystemBase {

  // Data members:

  // Time step ID
  int m_time_step;

  // Overflow counter and maximum number of steps
  int m_overflow_counter;
  int m_overflow_limit = std::numeric_limits<int>::max();

  // Analysis time
  Real m_time;

  // Stable time step size
  Real m_dt_scale_factor;
  Real m_stable_dt;

  // Global body forces
  Real m_bx = 0.0;
  Real m_by = 0.0;

  // Mass-proportional damping factor
  Real m_alpha = 0.0;

  // Contact stiffness parameter
  Real m_contact_stiffness = 0.0;

  // Initial velocity
  Real m_vx0 = 0.0;
  Real m_vy0 = 0.0;

  // Mesh object and data
  //Mesh m_mesh;

  // Boundary conditions
  //std::vector<BoundaryCondition> m_bcs;
  //std::vector<NodalForce> m_forces;

  // Node data
  int                     Nnodes; // The total number of nodes
  int             Ndofs_per_node; // The total number of degrees of freedom per node
  int                      Ndofs; // The total number of system degrees of freedom = Nnodes * Ndofs_per_node
  std::vector<bool>       fixity; // The fixity (fixed == true) assigned to each nodal degree of freedom
  std::vector<Real>            x; // The (primal) system initial position degrees of freedom
  std::vector<Real>           xt; // The (primal) system current position degrees of freedom
  std::vector<Dual<FixedU> >   u; // The (primal/dual) system displacement degrees of freedom
  std::vector<Dual<FixedV> >   v; // The (primal/dual) system velocity degrees of freedom
  std::vector<std::vector<FixedV> > v_overflow; // Additional storage for velocity overflow
  std::vector<Real>            m; // The (primal) system masses for each DoF
  std::vector<Real>            f; // The (primal) system forces for each DoF
  std::vector<Dual<Real> >     a; // The (primal) system accelerations for each DoF
  std::vector<Real>        alpha; // The system damping factors for each DoF
  std::vector<std::string> node_field_names;

  // Solid element data
  Element_T            m_element; // Solid element class
  int                     Nelems; // The total number of solid elements
  int            Nnodes_per_elem; // The total number of nodes per solid element
  std::vector<int>       connect; // The nodal connectivity array for all solid elements
  std::vector<Real>        state; // Solid element state variable data
  std::vector<std::string> element_field_names;

  // Truss element data
  Truss_T                m_truss; // Truss element class
  int                 Ntruss = 0; // The total number of truss elements
  std::vector<int> truss_connect; // The nodal connectivity array for all truss elements
  std::vector<Real>  truss_state; // Truss element state variable data
  std::vector<std::vector<Real> > truss_state_overflow; // Truss element state variable overflow data
  std::vector<std::string> truss_field_names;

  // Point mass data
  int              Npoints = 0; // The total number of point masses
  std::vector<int>   point_ids; // Node IDs of point massess
  std::vector<Real> point_mass; // Discrete mass associated with correponding nodal points
  std::vector<std::string> point_field_names;

  // Total system state data
  Real elastic_strain_energy;
  Real kinetic_energy;
  Real potential_energy;
  Real total_energy;
  std::vector<std::string> global_field_names;
  
  Integrator<Ratio>               m_integrator;           // (Bit-reversible) leapfrog time integrator
  std::vector<ContactInteraction> m_contact_interactions; // List of penalty-based contact interactions
    
  // ===================================================================== //
  
  // default constructor method for a system object
  System() : SystemBase() { }
    
  // ===================================================================== //

  // Procedure to initialize the system
  virtual void initialize(double* new_coordinates, double* new_velocities, bool* new_fixity, int new_Nnodes, int new_Ndofs_per_node,
		          int* new_connectivity, int new_Nelems, int new_Nnodes_per_elem,
		          Parameters& params) {
    
    std::cout << "| ====================== INITIALIZING ====================== |" << std::endl;

    std::cout << "Initializing state data..." << std::endl;

    // Initialize global parameters
    m_time_step = 0;
    m_time = 0.0;
    if (params.count("body_force_x") > 0) m_bx = params["body_force_x"];
    if (params.count("body_force_y") > 0) m_by = params["body_force_y"];
    if (params.count("initial_velocity_x") > 0) m_vx0 = params["initial_velocity_x"];
    if (params.count("initial_velocity_y") > 0) m_vy0 = params["initial_velocity_y"];
    if (params.count("dt_scale_factor") > 0) m_dt_scale_factor = params["dt_scale_factor"];
    if (params.count("mass_damping_factor") > 0) m_alpha = params["mass_damping_factor"];
    if (params.count("contact_stiffness") > 0) m_contact_stiffness = params["contact_stiffness"];
    if (params.count("overflow_limit") > 0) m_overflow_limit = int(params["overflow_limit"]);

    // Initialize the element/material object
    m_element = Element_T(params);

    // Initialize all mesh totals
    Nnodes          = new_Nnodes;
    Ndofs_per_node  = new_Ndofs_per_node;
    Ndofs           = Nnodes*Ndofs_per_node;
    Nelems          = new_Nelems;
    Nnodes_per_elem = new_Nnodes_per_elem;
    const int Nelem_dofs = Nelems*Nnodes_per_elem;
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    const int Nstate = Nelems*Nstate_vars_per_elem;
    if (Nelems > 0) element_field_names = m_element.m_model.get_field_names();

    // Initialize the dimensions of all arrays
    connect.resize(Nelem_dofs);
    fixity.resize(Ndofs);
    x.resize(Ndofs);
    xt.resize(Ndofs);
    u.resize(Ndofs,Dual<FixedU>(0.0,0.0));
    v.resize(Ndofs,Dual<FixedV>(0.0,0.0));
    m.resize(Ndofs);
    f.resize(Ndofs);
    a.resize(Ndofs,Dual<Real>(0.0,0.0));
    alpha.resize(Ndofs);
    state.resize(Nstate);

    // Populate ordered nodal field names
    node_field_names.push_back("displacement_X");
    node_field_names.push_back("displacement_Y");
    node_field_names.push_back("displacement_Z");
    node_field_names.push_back("velocity_X");
    node_field_names.push_back("velocity_Y");
    node_field_names.push_back("velocity_Z");
    node_field_names.push_back("force_X");
    node_field_names.push_back("force_Y");
    node_field_names.push_back("force_Z");
    node_field_names.push_back("dual_displacement_X");
    node_field_names.push_back("dual_displacement_Y");
    node_field_names.push_back("dual_displacement_Z");
    node_field_names.push_back("dual_velocity_X");
    node_field_names.push_back("dual_velocity_Y");
    node_field_names.push_back("dual_velocity_Z");

    // Initialize data for all DoFs
    for (int i=0; i<Ndofs; i++) {
      fixity[i] = new_fixity[i];
      x[i] = new_coordinates[i];
      xt[i] = x[i];
      u[i] = Dual<FixedU>(0.0,0.0);
      v[i] = Dual<FixedV>(new_velocities[i],0.0);
      m[i] = 0.0;
      f[i] = 0.0;
      a[i] = Dual<Real>(0.0,0.0);
      alpha[i] = 0.0;
    }
    
    // Assign constant initial velocity (if defined)
    for (int i = 0; i < Nnodes; i++) {
      if (m_vx0 != 0.0) v[2*i+0] = Dual<FixedV>(m_vx0,0.0);
      if (m_vy0 != 0.0) v[2*i+1] = Dual<FixedV>(m_vy0,0.0);
    }

    // Initialize dual DoFs
    for (int i=0; i<Ndofs; i++) {
      v[i].second = smallest_value(v[i].first);
    }

    // Initialize connectivity data for all elements
    for (int i=0; i<Nelem_dofs; i++) {
      connect[i] = new_connectivity[i];
    }

    // Initialize state variable data for all elements
    for (int e=0; e<Nelems; e++) {
      m_element.initialize(&state[Nstate_vars_per_elem*e]);
    }

    // Initialize system state data
    elastic_strain_energy = 0.0; global_field_names.push_back("elastic_strain_energy");
    kinetic_energy        = 0.0; global_field_names.push_back("kinetic_energy");
    potential_energy      = 0.0; global_field_names.push_back("potential_energy");
    total_energy          = 0.0; global_field_names.push_back("total_energy");
    
    std::cout << "| ========================================================== |" << std::endl;
    
  } // initialize()
  
  // ===================================================================== //

  // Procedure to initialize a new contact interaction
  virtual void initialize_contact(int* node_ids, int* segment_connectivity, int new_Nnodes, int new_Nsegments, Parameters& params) {
    
    // Define a new contact interaction, and initialize it
    m_contact_interactions.push_back(ContactInteraction());
    m_contact_interactions.back().initialize(node_ids,segment_connectivity,new_Nnodes,new_Nsegments,params);
    
  } // initialize_contact()
  
  // ===================================================================== //

  // Procedure to initialize truss elements
  virtual void initialize_truss_elements(int* new_truss_connectivity, int new_Ntruss, Parameters& params) {
    
    std::cout << "| =============== INITIALIZING TRUSS ELEMENTS ============== |" << std::endl;

    std::cout << "Initializing truss element state data..." << std::endl;

    // Initialize the truss element/material object
    m_truss = Truss_T(params);

    // Initialize truss totals
    Ntruss = new_Ntruss;
    const int Nnodes_per_truss = 2;
    const int Ntruss_dofs = Ntruss*Nnodes_per_truss;
    const int Nstate_vars_per_truss = m_truss.num_state_vars();
    const int Nstate = Ntruss*Nstate_vars_per_truss;
    if (Ntruss > 0) truss_field_names = m_truss.m_model.get_field_names();

    // Initialize the dimensions of all arrays
    truss_connect.resize(Ntruss_dofs);
    truss_state.resize(Nstate);
    truss_state_overflow.resize(Ntruss);

    // Initialize connectivity data for all truss elements
    for (int i=0; i<Ntruss_dofs; i++) {
      truss_connect[i] = new_truss_connectivity[i];
    }

    // Initialize state variable data for all truss elements
    for (int e=0; e<Ntruss; e++) {
      m_truss.initialize(&truss_state[Nstate_vars_per_truss*e]);
    }
    
    std::cout << "| ========================================================== |" << std::endl;
    
  } // initialize_truss_elements()
  
  // ===================================================================== //

  // Procedure to initialize new point masses
  virtual void initialize_point_mass(int* new_point_ids, double* new_point_mass, int new_Npoints, Parameters& params) {
    
    std::cout << "| ================ INITIALIZING POINT MASSES =============== |" << std::endl;

    std::cout << "Initializing point masses..." << std::endl;

    // Initialize point mass totals
    Npoints = new_Npoints;

    // Initialize the dimensions of all arrays
    point_ids.resize(Npoints);
    point_mass.resize(Npoints);

    // Initialize point mass data for all points
    for (int i=0; i<Npoints; i++) {
      point_ids[i]  = new_point_ids[i];
      point_mass[i] = new_point_mass[i];
    }
    if (Npoints > 0) point_field_names.push_back("mass");
    
    std::cout << "| ========================================================== |" << std::endl;
    
  } // initialize_contact()
  
  // ===================================================================== //

  // Procedure to initialize variable material stiffness properties
  virtual void initialize_variable_properties(double (*function_xy)(double,double)) {

    const int Nstate_vars_per_elem = m_element.num_state_vars();

    // Loop over all solid elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal positions for each element
      const int Ndofs_per_elem = Nnodes_per_elem*Ndofs_per_node;
      Real xe[8];
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
	}
      }

      // Initialize variable material properties for each element
      m_element.initialize_variable_properties(xe,&state[Nstate_vars_per_elem*e],function_xy);

    } // End loop over all solid elements
    
  } // initialize_variable_properties()
  
  // ===================================================================== //

  // Procedure to initialize the system state at the indicated time
  virtual void initialize_state(void) {

    // Set the initial time to zero
    m_time = 0.0;

    // Set the overflow counter
    m_overflow_counter = 0;
    
    // Update accelerations and damping factors for each DoF
    update_accelerations(0.0);

    // Update kinetic energy
    update_kinetic_energy();
    
  } // initialize_state()
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  virtual double update_state(Real dt) {

    const int Nstate_vars_per_truss = m_truss.num_state_vars();

    // Update the time step ID and conditionally load overflow dual velocities
    if (dt < 0.0) {
      m_time_step--;
      if (m_overflow_counter == 0) {
	m_overflow_counter = m_overflow_limit;
	std::vector<FixedV>& last_v_overflow = v_overflow.back();
	for (int i=0; i<Ndofs; i++) {
	  v[i].second = last_v_overflow[i];
	}
	v_overflow.pop_back();
	std::cout << "Loading velocity overflow: count = " << v_overflow.size() << std::endl;
      }
      m_overflow_counter--;

      // Conditionally load material history parameters from memory
      for (int e=0; e<Ntruss; e++) {
	if (truss_state_overflow[e].size() > 0) {
	  if (m_truss.load_state(&truss_state[Nstate_vars_per_truss*e],&truss_state_overflow[e].back())) {
	    truss_state_overflow[e].pop_back();
	    std::cout << "Loading state overflow: count = " << truss_state_overflow[e].size() << std::endl;
	  }
	}
      }
      
    }

    // Update velocities to the half-step
    m_integrator.first_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),Ndofs);

    // Update displacement to the next whole-step
    m_integrator.whole_step_displacement_update(dt,v.data(),u.data(),Ndofs);
    
    // Update masses, residual forces, and accelerations at the whole-step
    update_accelerations(dt);
    
    // Update velocities to the whole-step
    m_integrator.second_half_step_velocity_update(dt,v.data(),a.data(),alpha.data(),Ndofs);

    // Update kinetic energy
    update_kinetic_energy();

    // Update time step ID and conditionally store overflow dual velocities
    if (dt > 0.0) {
      m_time_step++;
      m_overflow_counter++;
      if (m_overflow_counter == m_overflow_limit) {
	m_overflow_counter = 0;
	std::vector<FixedV> new_v_overflow(Ndofs,FixedV(0.0));
	for (int i=0; i<Ndofs; i++) {
	  //if ((v[i].first == FixedV(0.0)) and !fixity[i]) exit(1);
	  new_v_overflow[i] = v[i].second;
	  v[i].second = smallest_value(v[i].second);
	}
	v_overflow.push_back(new_v_overflow);
	std::cout << "Storing velocity overflow: count = " << v_overflow.size() << std::endl;
      }

      // Conditionally store material history parameters in memory
      for (int e=0; e<Ntruss; e++) {
	Real new_state_overflow;
	if (m_truss.store_state(&truss_state[Nstate_vars_per_truss*e],&new_state_overflow)) {
	  truss_state_overflow[e].push_back(new_state_overflow);
	  std::cout << "Storing state overflow: count = " << truss_state_overflow[e].size() << std::endl;
	}
      }
      
    }
    
    std::cout << "Time step: " << m_time_step << " at time: " << m_time << std::endl;

    // Return the updated analysis time
    return m_time;
    
  } // update_state()
  
  // ===================================================================== //

  // Procedure to get the current system state data
  virtual double get_field_data(double *ux, double *uy,
	 	                double *vx, double *vy,
		                double *fx, double *fy,
			        double *dual_ux, double *dual_uy,
	 	                double *dual_vx, double *dual_vy,
		                double *sxx, double *syy, double *sxy,
				double *pressure, double *stiffness_scaling_factor,
				double *system_state, double *eqps, bool *is_dead) {

    // Copy nodal state data
    for (int i=0; i<Nnodes; i++) {
      ux[i]      = u[2*i+0].first;
      uy[i]      = u[2*i+1].first;
      vx[i]      = v[2*i+0].first;
      vy[i]      = v[2*i+1].first;
      fx[i]      = f[2*i+0];
      fy[i]      = f[2*i+1];
      dual_ux[i] = u[2*i+0].second;
      dual_uy[i] = u[2*i+1].second;
      dual_vx[i] = v[2*i+0].second;
      dual_vy[i] = v[2*i+1].second;
    }

    // Copy element state data
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    for (int e=0; e<Nelems; e++) {
      sxx[e]                      = state[Nstate_vars_per_elem*e+0];
      syy[e]                      = state[Nstate_vars_per_elem*e+1];
      sxy[e]                      = state[Nstate_vars_per_elem*e+2];
      pressure[e]                 = state[Nstate_vars_per_elem*e+3];
      stiffness_scaling_factor[e] = state[Nstate_vars_per_elem*e+4];
    }

    // Copy system state data
    system_state[0] = elastic_strain_energy;
    system_state[1] = kinetic_energy;
    system_state[2] = potential_energy;
    system_state[3] = total_energy;
    
    // Copy truss state data
    const int Nstate_vars_per_truss = m_truss.num_state_vars();
    for (int e=0; e<Ntruss; e++) {
      eqps[e]    = m_truss.get_state_variable(&truss_state[Nstate_vars_per_truss*e],"eqps");
      is_dead[e] = m_truss.is_dead(&truss_state[Nstate_vars_per_truss*e]);
    }

    // Return the current analysis time
    return m_time;
    
  } // get_field_data()
  
  // ===================================================================== //
  
  // Request the number of spatial dimensions
  virtual int get_num_dim(void) { return 2; }
  
  // ===================================================================== //

  // Request the number of entities of the specified type
  virtual int get_num_entities(std::string entity_type) {
    if      (entity_type == "global")  { return 1;       }
    else if (entity_type == "node")    { return Nnodes;  }
    else if (entity_type == "element") { return Nelems;  }
    else if (entity_type == "truss")   { return Ntruss;  }
    else if (entity_type == "point")   { return Npoints; }
    else                               { return 0;       }
  } // get_num_entitites()
  
  // ===================================================================== //

  // Request the number of fields defined for entities of the specified type
  virtual int get_num_fields(std::string entity_type) {
    if      (entity_type == "global")  { return global_field_names.size();  }
    else if (entity_type == "node")    { return node_field_names.size();    }
    else if (entity_type == "element") { return element_field_names.size(); }
    else if (entity_type == "truss")   { return truss_field_names.size();   }
    else if (entity_type == "point")   { return point_field_names.size();   }
    else                               { return 0;                          }
  } // get_num_fields()
  
  // ===================================================================== //

  // Request the name of the indicated field ID for entities of the specified type
  virtual const char* get_field_name(std::string entity_type, int field_id) {
    if      (entity_type == "global")  { return global_field_names[field_id].c_str();  }
    else if (entity_type == "node")    { return node_field_names[field_id].c_str();    }
    else if (entity_type == "element") { return element_field_names[field_id].c_str(); }
    else if (entity_type == "truss")   { return truss_field_names[field_id].c_str();   }
    else if (entity_type == "point")   { return point_field_names[field_id].c_str();   }
    else                               { return nullptr;                               }
  } // get_field_name()
  
  // ===================================================================== //
  
  // Request the coordinates of all nodes
  virtual void get_node_coords(double* coords, bool deformed) {
    if (deformed) { // return the currently deformed coordinates of all nodes
      for (int i=0; i<Ndofs; i++) { coords[i] = xt[i]; }
    } else {        // return the initial (undeformed) coordinates of all nodes
      for (int i=0; i<Ndofs; i++) { coords[i] = x[i];  }
    }
  } // get_node_coords()
  
  // ===================================================================== //
  
  // Request the connectivity data for all entities of the specified type
  virtual void get_connectivity(std::string entity_type, int* connectivity) {
    if        (entity_type == "element") {
      for (int i=0; i<Nelems; i++) {
	for (int j=0; j<Nnodes_per_elem; j++) {
	  connectivity[Nnodes_per_elem*i+j] = connect[Nnodes_per_elem*i+j];
	}
      }
    } else if (entity_type == "truss")   {
      for (int i=0; i<Ntruss; i++) {
	for (int j=0; j<2; j++) {
	  connectivity[2*i+j] = truss_connect[2*i+j];
	}
      }
    } else if (entity_type == "point")   {
      for (int i=0; i<Npoints; i++) {
	connectivity[i] = point_ids[i];
      }
    }
  } // get_connectivity()
  
  // ===================================================================== //

  // Request data defining all fields for all entities of the specified type
  virtual void get_fields(std::string entity_type, double* field_data) {
    if        (entity_type == "global")  {
      field_data[0] = elastic_strain_energy;
      field_data[1] = kinetic_energy;
      field_data[2] = potential_energy;
      field_data[3] = total_energy;
    } else if (entity_type == "node")    {
      const int Nstate = get_num_fields(entity_type);
      for (int i=0; i<Nnodes; i++) {
	field_data[Nstate*i+ 0] = u[2*i+0].first;
	field_data[Nstate*i+ 1] = u[2*i+1].first;
	field_data[Nstate*i+ 2] = 0.0;
	field_data[Nstate*i+ 3] = v[2*i+0].first;
	field_data[Nstate*i+ 4] = v[2*i+1].first;
	field_data[Nstate*i+ 5] = 0.0;
	field_data[Nstate*i+ 6] = f[2*i+0];
	field_data[Nstate*i+ 7] = f[2*i+1];
	field_data[Nstate*i+ 8] = 0.0;
	field_data[Nstate*i+ 9] = u[2*i+0].second;
	field_data[Nstate*i+10] = u[2*i+1].second;
	field_data[Nstate*i+11] = 0.0;
	field_data[Nstate*i+12] = v[2*i+0].second;
	field_data[Nstate*i+13] = v[2*i+1].second;
	field_data[Nstate*i+14] = 0.0;
      }
    } else if (entity_type == "element") {
      const int Nstate = m_element.num_state_vars();
      const int Nmat_state = get_num_fields(entity_type);
      for (int e=0; e<Nelems; e++) {
        m_element.m_model.get_fields(&state[Nstate*e],&field_data[Nmat_state*e]);
      }
    } else if (entity_type == "truss")   {
      const int Nstate = m_truss.num_state_vars();
      const int Nmat_state = m_truss.m_model.num_state_vars();
      for (int e=0; e<Ntruss; e++) {
        m_truss.m_model.get_fields(&truss_state[Nstate*e],&field_data[Nmat_state*e]);
      }
    } else if (entity_type == "point")   {
      const int Nstate = get_num_fields(entity_type);
      for (int i=0; i<Npoints; i++) {
	field_data[Nstate*i+0] = point_mass[i];
      }
    }
  } // get_fields()
  
  // ===================================================================== //

  // Request the current analysis time
  virtual double get_time(void) { return m_time; }
  
  // ===================================================================== //
private:
  // ===================================================================== //

  // Procedure to update the accelerations
  void update_accelerations(Real dt) {

    // Update the current analysis time
    m_time = m_time + dt;

    // Zero-initialize forces and masses
    std::fill(m.begin(), m.end(), 0.0);
    std::fill(f.begin(), f.end(), 0.0);

    // Zero-initialize the total system elastic strain energy and potential energy
    elastic_strain_energy = 0.0;
    potential_energy      = 0.0;

    // Update the current deformed nodal coordinates
    for (int i = 0; i < Nnodes; i++) {
      xt[2*i+0] = x[2*i+0] + Real(u[2*i+0].first);
      xt[2*i+1] = x[2*i+1] + Real(u[2*i+1].first);
    }
    
    const int Nstate_vars_per_elem = m_element.num_state_vars();

    // Loop over all solid elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal displacements for each element (cast as Reals)
      const int Ndofs_per_elem = Nnodes_per_elem*Ndofs_per_node;
      Real xe[8];
      Real ue[8];
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
	  ue[Ndofs_per_node*j+i] = u[Ndofs_per_node*jnode_id+i].first;
	}
      }

      // Compute internal force in the current element
      Real me[8] = { 0.0 };
      Real fe[8] = { 0.0 };
      Real Ee = 0.0;
      m_element.update(xe,ue,me,fe,Ee,&state[Nstate_vars_per_elem*e],dt);

      // Scatter mass and forces to the nodes
      // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
      for (int j=0; j<Nnodes_per_elem; j++) {
	const int jnode_id = connect[Nnodes_per_elem*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  m[Ndofs_per_node*jnode_id+i] += me[Ndofs_per_node*j+i];
	  f[Ndofs_per_node*jnode_id+i] -= fe[Ndofs_per_node*j+i];
	}
      }

      // Sum contribution to the total elastic strain energy
      // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
      elastic_strain_energy += Ee;

    } // End loop over all solid elements
    
    const int Nnodes_per_truss = 2;
    const int Nstate_vars_per_truss = m_truss.num_state_vars();

    // Loop over all truss elements
    for (int e=0; e<Ntruss; e++) {

      // Copy (primal) local nodal displacements for each truss element (cast as Reals)
      const int Ndofs_per_truss = Nnodes_per_truss*Ndofs_per_node;
      Real xe[4];
      Real ue[4];
      for (int j=0; j<Nnodes_per_truss; j++) {
	const int jnode_id = truss_connect[Nnodes_per_truss*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
	  ue[Ndofs_per_node*j+i] = u[Ndofs_per_node*jnode_id+i].first;
	}
      }

      // Compute internal force in the current truss element
      Real me[4] = { 0.0 };
      Real fe[4] = { 0.0 };
      Real Ee = 0.0;
      m_truss.update(xe,ue,me,fe,Ee,&truss_state[Nstate_vars_per_truss*e],dt);

      // Scatter mass and forces to the nodes
      // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
      for (int j=0; j<Nnodes_per_truss; j++) {
	const int jnode_id = truss_connect[Nnodes_per_truss*e+j];
	for (int i=0; i<Ndofs_per_node; i++) {
	  m[Ndofs_per_node*jnode_id+i] += me[Ndofs_per_node*j+i];
	  f[Ndofs_per_node*jnode_id+i] -= fe[Ndofs_per_node*j+i];
	}
      }

      // Sum contribution to the total elastic strain energy
      // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
      elastic_strain_energy += Ee;

    } // End loop over all truss elements

    // Loop over all point masses
    for (int i=0; i<Npoints; i++) {
      // sum point mass to correponding nodal DoFs
      m[2*point_ids[i]+0] += point_mass[i];
      m[2*point_ids[i]+1] += point_mass[i];
    } // End loop over all point masses
    
    // Sum nodal forces from external nodal point loads (gravity, etc.)
    for (int i = 0; i < Nnodes; i++) {
      // Global body forces
      f[2*i+0] += m_bx*m[2*i+0];
      f[2*i+1] += m_by*m[2*i+1];

      // Sum contributions to the total potential energy
      potential_energy -= (m_bx*m[2*i+0]*xt[2*i+0] + m_by*m[2*i+1]*xt[2*i+1]);
    }

    // Sum nodal forces due to contact interactions
    for (auto& contact : m_contact_interactions) {
      // Update contact forces
      Real Ec = 0.0;
      contact.update_contact_forces(xt.data(),f.data(),Ec,dt);

      // Sum contribution to the total elastic strain energy
      elastic_strain_energy += Ec;
    }

    // Contact interaction with rigid wall
    for (int i = 0; i < Nnodes; i++) {
      Real yt = x[2*i+1] + Real(u[2*i+1].first);
      if (yt < 0.0) {
	Real du = - yt;
	f[2*i+1] += m_contact_stiffness*du;

        // Sum contribution to the total elastic strain energy
	elastic_strain_energy += 0.5*m_contact_stiffness*du*du;
      }
    }

    // Compute updated accelerations and damping factors for each DoF:
    for (int i=0; i<Ndofs; i++) {

      // Divide nodal forces by nodal masses to obtain the (predictor) nodal accelerations
      a[i].first = f[i]/m[i];

      // Determine fictitious dual nodal accelerations (directly proportional to the dual displacements)
      a[i].second = 0.0*Real(u[i].second);

      // Determine the (mass-proportional) damping factor "alpha" in terms of a target frequency and damping ratio
      //Real frequency = 1.0;
      //Real damping_ratio = 0.1;
      //alpha[i] = 2.0*frequency*damping_ratio;
      alpha[i] = m_alpha;

      // Impose fixed boundary conditions
      if (fixity[i]) {
	a[i].first  = 0.0;
	a[i].second = 0.0;
	v[i].first  = 0.0;
	v[i].second = 0.0;
      }

    } // End loop over all DoFs

    // Report the largest dual variable value (indicates an integer overflow ...)
    //Integer max_val;
    //for (int i = 0; i < Ndofs; i++) {
    //  max_val = std::max(max_val,std::abs(v[i].second.mantissa));
    //}
    //if (Real(max_val) > 0.99*std::numeric_limits<Integer>::max()) std::cout << "max val = " << max_val << " > " << 0.99*std::numeric_limits<Integer>::max() << std::endl;
    
  } // update_accelerations()
  
  // ===================================================================== //

  // Determine the updated system kinetic energy
  void update_kinetic_energy() {

    // Zero-initialize the total system kinetic energy
    kinetic_energy = 0.0;

    // Loop over all DoFs and sum contributions to the total kinetic energy
    for (int i=0; i<Ndofs; i++) {
      Real vi = v[i].first;
      kinetic_energy += 0.5*m[i]*vi*vi;
    } // End loop over all DoFs

    // Update the total system energy
    total_energy = elastic_strain_energy + kinetic_energy + potential_energy;
    
  } // update_kinetic_energy()
  
  // ===================================================================== //
  
}; // System

#endif // SYSTEM_H

