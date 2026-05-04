#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.h"
#include "BoundaryCondition.h"
#include "ContactInteraction.h"
#include "Dual.h"
#include "Integrator.h"
#include "MaterialAdjoint.h"
#include "PassPhase.h"
#include "Parameters.h"
#include <limits>
#include <vector>
#include <algorithm> // For std::fill
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <math.h>
#include <stdlib.h> // exit

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

  // Procedure to define a time-varying displacement boundary condition
  virtual void define_displacement_bc(int* node_ids, int num_nodes, int component, TimeFunction function) = 0;
  
  // ===================================================================== //

  // Procedure to initialize variable material stiffness properties
  virtual void initialize_variable_properties(double (*function_xy)(double,double)) = 0;

  // Procedure to initialize variable material relaxation-time properties
  virtual void initialize_variable_relaxation_time(double (*function_xy)(double,double)) = 0;
  
  // ===================================================================== //

  // Procedure to initialize the system state at the indicated time
  virtual void initialize_state(void) = 0;
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  virtual double update_state(Real dt, PassPhase phase) = 0;
  
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

  // Clear global adjoint state (u*, v*, kick seeds and global adjoint accumulators)
  virtual void clear_adjoint_state(void) = 0;

  // ===================================================================== //

  // Add nodal velocity adjoint seeds (seed_xy stores [seed_x, seed_y] per node)
  virtual void add_nodal_velocity_adjoint_seed(const int* node_ids, const double* seed_xy, int num_nodes) = 0;

  // ===================================================================== //

  // Add nodal displacement adjoint seeds (seed_xy stores [seed_x, seed_y] per node)
  virtual void add_nodal_displacement_adjoint_seed(const int* node_ids, const double* seed_xy, int num_nodes) = 0;

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
  std::vector<Real> m_dt_history;
  std::vector<Real> m_time_history;

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
  std::vector<DisplacementBoundaryCondition> m_displacement_bcs; // Time-varying displacement boundary conditions
  std::vector<bool> m_has_time_bc; // Flags indicating DoFs with time-varying BCs

  // Solid element data
  Element_T            m_element; // Solid element class
  int                     Nelems; // The total number of solid elements
  int            Nnodes_per_elem; // The total number of nodes per solid element
  std::vector<int>       connect; // The nodal connectivity array for all solid elements
  std::vector<Real>        state; // Solid element state variable data
  std::vector<std::string> element_field_names;
  std::vector<std::vector<Real> > element_state_overflow; // Element integration-point overflow storage

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
  std::vector<Real> u_adjoint;
  std::vector<Real> v_adjoint;
  std::vector<Real> kick_force_seed;
  Real adjoint_mass_damping_gradient = 0.0;
  Real adjoint_contact_stiffness_gradient = 0.0;
  std::vector<std::string> adjoint_param_names;
  std::vector<Real> adjoint_param_gradients;

  // Optional adjoint diagnostics for stiffness-scaling gradient tracing.
  bool m_adjoint_debug_dump_enable = false;
  Real m_adjoint_debug_dump_threshold = 0.0;
  int m_adjoint_debug_dump_max_rows = 200000;
  int m_adjoint_debug_dump_stride = 1;
  int m_adjoint_debug_dump_rows = 0;
  LongInteger m_adjoint_debug_dump_run_id = 0;
  std::ofstream m_adjoint_debug_dump_stream;
  int m_element_stiffness_param_id = -1;

  // Optional overflow-risk warnings for dual/ancillary variables.
  bool m_dual_overflow_warn_enable = false;
  Real m_dual_overflow_warn_fraction = 0.95;
  int m_dual_overflow_warn_limit = 20;
  int m_dual_overflow_warn_count = 0;
  
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
    m_dt_history.clear();
    m_time_history.clear();
    node_field_names.clear();
    element_field_names.clear();
    truss_field_names.clear();
    point_field_names.clear();
    global_field_names.clear();
    adjoint_param_names.clear();
    adjoint_param_gradients.clear();
    adjoint_mass_damping_gradient = 0.0;
    adjoint_contact_stiffness_gradient = 0.0;
    if (params.count("body_force_x") > 0) m_bx = params["body_force_x"];
    if (params.count("body_force_y") > 0) m_by = params["body_force_y"];
    if (params.count("initial_velocity_x") > 0) m_vx0 = params["initial_velocity_x"];
    if (params.count("initial_velocity_y") > 0) m_vy0 = params["initial_velocity_y"];
    if (params.count("dt_scale_factor") > 0) m_dt_scale_factor = params["dt_scale_factor"];
    if (params.count("mass_damping_factor") > 0) m_alpha = params["mass_damping_factor"];
    if (params.count("contact_stiffness") > 0) m_contact_stiffness = params["contact_stiffness"];
    if (params.count("overflow_limit") > 0) m_overflow_limit = int(params["overflow_limit"]);
    if (params.count("adjoint_debug_dump_enable") > 0) m_adjoint_debug_dump_enable = (params["adjoint_debug_dump_enable"] != 0.0);
    if (params.count("adjoint_debug_dump_threshold") > 0) m_adjoint_debug_dump_threshold = params["adjoint_debug_dump_threshold"];
    if (params.count("adjoint_debug_dump_max_rows") > 0) m_adjoint_debug_dump_max_rows = std::max(0,int(params["adjoint_debug_dump_max_rows"]));
    if (params.count("adjoint_debug_dump_stride") > 0) m_adjoint_debug_dump_stride = std::max(1,int(params["adjoint_debug_dump_stride"]));
    if (params.count("dual_overflow_warn_enable") > 0) m_dual_overflow_warn_enable = (params["dual_overflow_warn_enable"] != 0.0);
    if (params.count("dual_overflow_warn_fraction") > 0) {
      m_dual_overflow_warn_fraction = std::max(0.0,std::min(1.0,params["dual_overflow_warn_fraction"]));
    }
    if (params.count("dual_overflow_warn_limit") > 0) m_dual_overflow_warn_limit = std::max(0,int(params["dual_overflow_warn_limit"]));
    m_dual_overflow_warn_count = 0;
    m_adjoint_debug_dump_rows = 0;
    m_adjoint_debug_dump_run_id += 1;

    // Initialize the element/material object
    m_element = Element_T(params);
    m_element_stiffness_param_id = -1;
    for (int p=0; p<m_element.adjoint_num_params(); p++) {
      if (std::string(m_element.adjoint_param_name(p)) == "stiffness_scaling_factor") {
        m_element_stiffness_param_id = p;
        break;
      }
    }
    initialize_adjoint_debug_dump_stream();

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

    // Reset displacement boundary condition data
    m_displacement_bcs.clear();
    m_has_time_bc.resize(Ndofs);

    // assign() → // initialize all DOFs (assign is cleaner than resize+manual loop)
    // Initialize the dimensions of all arrays
    connect.resize(Nelem_dofs);
    fixity.resize(Ndofs);
    x.resize(Ndofs);
    xt.resize(Ndofs);
    u.resize(Ndofs,Dual<FixedU>(0.0,0.0));
    v.resize(Ndofs,Dual<FixedV>(0.0,0.0));
    u_adjoint.resize(Ndofs,0.0);
    v_adjoint.resize(Ndofs,0.0);
    kick_force_seed.resize(Ndofs,0.0);
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
    node_field_names.push_back("adjoint_displacement_X");
    node_field_names.push_back("adjoint_displacement_Y");
    node_field_names.push_back("adjoint_displacement_Z");
    node_field_names.push_back("adjoint_velocity_X");
    node_field_names.push_back("adjoint_velocity_Y");
    node_field_names.push_back("adjoint_velocity_Z");

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
      m_has_time_bc[i] = false;
    }
    
    // Assign constant initial velocity (if defined)
    for (int i = 0; i < Nnodes; i++) {
      if (m_vx0 != 0.0) v[2*i+0] = Dual<FixedV>(m_vx0,0.0);
      if (m_vy0 != 0.0) v[2*i+1] = Dual<FixedV>(m_vy0,0.0);
    }

    // Initialize dual DoFs
    for (int i=0; i<Ndofs; i++) {
      v[i].second = FixedV(0.0);
    }

    // Initialize connectivity data for all elements
    for (int i=0; i<Nelem_dofs; i++) {
      connect[i] = new_connectivity[i];
    }

    // Initialize state variable data for all elements
    for (int e=0; e<Nelems; e++) {
      m_element.initialize(&state[Nstate_vars_per_elem*e]);
    }

    // Initialize element overflow storage
    element_state_overflow.resize(Nelems);

    // Initialize system state data
    elastic_strain_energy = 0.0; global_field_names.push_back("elastic_strain_energy");
    kinetic_energy        = 0.0; global_field_names.push_back("kinetic_energy");
    potential_energy      = 0.0; global_field_names.push_back("potential_energy");
    total_energy          = 0.0; global_field_names.push_back("total_energy");
    register_adjoint_params_from_materials();
    
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

    register_adjoint_params_from_materials();
    
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

  // Procedure to define a time-varying displacement boundary condition
  virtual void define_displacement_bc(int* node_ids, int num_nodes, int component, TimeFunction function) {
    
    if (function == nullptr) {
      std::cerr << "Null time function provided for displacement boundary condition." << std::endl;
      return;
    }

    if (Ndofs == 0) {
      std::cerr << "Geometry must be defined before adding displacement boundary conditions." << std::endl;
      return;
    }

    if ((component < 0) || (component >= Ndofs_per_node)) {
      std::cerr << "Invalid component (" << component << ") for displacement boundary condition." << std::endl;
      return;
    }

    DisplacementBoundaryCondition bc;
    bc.component = component;
    bc.function = function;
    bc.nodes.node_ids.clear();
    bc.last_values.assign(num_nodes, 0.0);
    bc.prescribed_velocities.assign(num_nodes, 0.0);
    
    // Loop over all nodes that will have the displacement BC applied
    for (int i=0; i<num_nodes; ++i) {
      int node_id = node_ids[i];
      if ((node_id < 0) || (node_id >= Nnodes)) {
	std::cerr << "Ignoring invalid node id (" << node_id << ") in displacement boundary condition." << std::endl;
	continue;
      }

      int dof = Ndofs_per_node*node_id + component;
      if ((dof < 0) || (dof >= Ndofs)) {
	std::cerr << "Ignoring invalid DoF index (" << dof << ") in displacement boundary condition." << std::endl;
	continue;
      }

      bc.nodes.node_ids.push_back(node_id);
      int constrained_dof = bc.nodes.node_ids.size() - 1;
      
      Real px = x[Ndofs_per_node*node_id + 0];
      Real py = 0.0;
      if (Ndofs_per_node > 1) py = x[Ndofs_per_node*node_id + 1];

      Real initial_value = function(m_time, px, py);
      bc.last_values[constrained_dof] = initial_value;
      bc.prescribed_velocities[constrained_dof] = 0.0;

      fixity[dof] = true;
      m_has_time_bc[dof] = true;
      u[dof] = Dual<FixedU>(initial_value, 0.0);
      v[dof] = Dual<FixedV>(0.0, 0.0);
    }

    if (!bc.nodes.node_ids.empty()) {
      int constrained_dof = bc.nodes.node_ids.size();
      bc.last_values.resize(constrained_dof);
      bc.prescribed_velocities.resize(constrained_dof);
      m_displacement_bcs.push_back(bc);
    }
    
  } // define_displacement_bc()
  
  // ===================================================================== //

  // Procedure to initialize variable material stiffness properties
  virtual void initialize_variable_properties(double (*function_xy)(double,double)) {

    const int Nstate_vars_per_elem = m_element.num_state_vars();

    // Loop over all solid elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal positions for each element
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

  // Procedure to initialize variable material relaxation-time properties
  virtual void initialize_variable_relaxation_time(double (*function_xy)(double,double)) {

    const int Nstate_vars_per_elem = m_element.num_state_vars();

    // Loop over all solid elements
    for (int e=0; e<Nelems; e++) {

      // Copy (primal) local nodal positions for each element
      Real xe[8];
      for (int j=0; j<Nnodes_per_elem; j++) {
        const int jnode_id = connect[Nnodes_per_elem*e+j];
        for (int i=0; i<Ndofs_per_node; i++) {
          xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
        }
      }

      // Initialize variable material relaxation-time values for each element
      m_element.initialize_variable_relaxation_time(xe,&state[Nstate_vars_per_elem*e],function_xy);

    } // End loop over all solid elements

  } // initialize_variable_relaxation_time()
  
  // ===================================================================== //

  // Procedure to initialize the system state at the indicated time
  virtual void initialize_state(void) {

    // Set the initial time to zero
    m_time = 0.0;

    // Set the overflow counter
    m_overflow_counter = 0;
    v_overflow.clear();
    for (auto& overflow_e : truss_state_overflow) {
      overflow_e.clear();
    }
    for (auto& overflow_e : element_state_overflow) {
      overflow_e.clear();
    }
    m_dt_history.clear();
    m_time_history.clear();
    std::fill(u_adjoint.begin(),u_adjoint.end(),0.0);
    std::fill(v_adjoint.begin(),v_adjoint.end(),0.0);
    std::fill(kick_force_seed.begin(),kick_force_seed.end(),0.0);
    adjoint_mass_damping_gradient = 0.0;
    adjoint_contact_stiffness_gradient = 0.0;

    // Apply any prescribed displacement boundary conditions at the initial time
    apply_displacement_bcs(m_time,0.0);
    enforce_prescribed_velocities();
    
    // Update accelerations and damping factors for each DoF
    update_accelerations(0.0,PassPhase::Forward,0.0);

    // Update kinetic energy
    update_kinetic_energy();
    update_global_param_gradients();
    
  } // initialize_state()
  
  // ===================================================================== //

  // Procedure to update the system state for a given time step
  virtual double update_state(Real dt, PassPhase phase) {

    if (dt <= 0.0) {
      std::cout << "System::update_state requires dt > 0.0; received dt = " << dt << std::endl;
      exit(EXIT_FAILURE);
    }

    const int Nstate_vars_per_truss = m_truss.num_state_vars();
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    const bool reverse_phase = is_reverse_phase(phase);

    Real dt_step = dt;
    Real signed_dt = dt;

    // Update the time step ID and conditionally load overflow dual velocities
    if (reverse_phase) {
      if (m_dt_history.empty()) {
        std::cout << "System::update_state reverse requested with empty time-step history" << std::endl;
        exit(EXIT_FAILURE);
      }

      m_time_step--;
      dt_step = m_dt_history.back();
      m_dt_history.pop_back();
      signed_dt = -dt_step;

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
      // Load truss overflow
      int Ntruss_overflow = 0;
      for (int e=0; e<Ntruss; e++) {
        int old_state_overflow_size = truss_state_overflow[e].size();
        m_truss.load_state(&truss_state[Nstate_vars_per_truss*e],truss_state_overflow[e]);
        if (truss_state_overflow[e].size() < old_state_overflow_size) { Ntruss_overflow++; }
      }
      if (Ntruss_overflow > 0) { std::cout << "Loaded " << Ntruss_overflow << " truss overflow states" << std::endl; }

      // Conditionally load material history parameters from memory
      // Load element overflow
      int Nelem_overflow = 0;
      for (int e=0; e<Nelems; e++) {
        int old_state_overflow_size = element_state_overflow[e].size();
        m_element.load_state(&state[Nstate_vars_per_elem*e],element_state_overflow[e]);
        if (element_state_overflow[e].size() < old_state_overflow_size) { Nelem_overflow++; }
      }
      if (Nelem_overflow > 0) { std::cout << "Loaded " << Nelem_overflow << " element overflow states" << std::endl; }
    }

    if (phase == PassPhase::BackwardAdjoint) {
      enforce_backward_adjoint_guardrails();
      apply_second_kick_damping_pullback(dt_step);
      // Second kick adjoint at the current reconstructed state (n+1).
      apply_global_kick_adjoint(dt_step,false);
    }

    // Update velocities to the half-step
    m_integrator.first_half_step_velocity_update(signed_dt,v.data(),a.data(),alpha.data(),Ndofs);

    // Update displacement to the next whole-step
    m_integrator.whole_step_displacement_update(signed_dt,v.data(),u.data(),Ndofs);

    // Storing/restoring time history
    if (signed_dt > 0.0) {
      m_time_history.push_back(m_time);
      m_time = m_time + signed_dt;
    } else if (signed_dt < 0.0) {
      m_time = m_time_history.back();
      m_time_history.pop_back();
    }

    // Enforce prescribed displacement boundary conditions at the new analysis time
    apply_displacement_bcs(m_time,signed_dt);

    if (phase == PassPhase::BackwardAdjoint) {
      // Drift transpose contribution (u* -> v*) after reverse drift.
      for (int i=0; i<Ndofs; i++) {
        v_adjoint[i] += u_adjoint[i]*dt_step;
      }
    }

    // Update masses, residual forces, and accelerations at the whole-step
    update_accelerations(signed_dt,phase,dt_step);
    // enforce_prescribed_velocities();
    // Update velocities to the whole-step
    m_integrator.second_half_step_velocity_update(signed_dt,v.data(),a.data(),alpha.data(),Ndofs);

    if (phase == PassPhase::BackwardAdjoint) {
      apply_first_kick_damping_pullback(dt_step);
    }

    maybe_warn_dual_velocity_overflow("post_second_half_step");
    enforce_prescribed_velocities();

    // Update kinetic energy
    update_kinetic_energy();

    // Update time step ID and conditionally store overflow dual velocities
    if (!reverse_phase) {
      m_dt_history.push_back(dt_step);
      m_time_step++;
      m_overflow_counter++;
      if (m_overflow_counter == m_overflow_limit) {
        m_overflow_counter = 0;
        std::vector<FixedV> new_v_overflow(Ndofs,FixedV(0.0));
        for (int i=0; i<Ndofs; i++) {
          //if ((v[i].first == FixedV(0.0)) and !fixity[i]) exit(1);
          new_v_overflow[i] = v[i].second;
          v[i].second = FixedV(0.0);
        }
        v_overflow.push_back(new_v_overflow);
        std::cout << "Storing velocity overflow: count = " << v_overflow.size() << std::endl;
      }

      // Conditionally store material history parameters in memory
      // Store truss overflow
      int Ntruss_overflow = 0;
      for (int e=0; e<Ntruss; e++) {
        int old_state_overflow_size = truss_state_overflow[e].size();
        m_truss.store_state(&truss_state[Nstate_vars_per_truss*e],truss_state_overflow[e]);
        if (truss_state_overflow[e].size() > old_state_overflow_size) { Ntruss_overflow++; }
      }
      if (Ntruss_overflow > 0) { std::cout << "Stored " << Ntruss_overflow << " truss overflow states" << std::endl; }

      // Conditionally store material history parameters in memory
      // Store element overflow
      int Nelem_overflow = 0;
      for (int e=0; e<Nelems; e++) {
        int old_state_overflow_size = element_state_overflow[e].size();
        m_element.store_state(&state[Nstate_vars_per_elem*e],element_state_overflow[e]);
        if (element_state_overflow[e].size() > old_state_overflow_size) { Nelem_overflow++; }
      }
      if (Nelem_overflow > 0) { std::cout << "Stored " << Nelem_overflow << " element overflow states" << std::endl; }
    }
    
    std::cout << "Time step: " << m_time_step << " at time: " << m_time << std::endl;

    update_global_param_gradients();

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
      for (int i=0; i<int(adjoint_param_gradients.size()); i++) {
        field_data[4+i] = adjoint_param_gradients[i];
      }
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
		field_data[Nstate*i+15] = u_adjoint[2*i+0];
		field_data[Nstate*i+16] = u_adjoint[2*i+1];
		field_data[Nstate*i+17] = 0.0;
		field_data[Nstate*i+18] = v_adjoint[2*i+0];
		field_data[Nstate*i+19] = v_adjoint[2*i+1];
		field_data[Nstate*i+20] = 0.0;
	      }
    } else if (entity_type == "element") {
      const int Nstate = m_element.num_state_vars();
      const int Nmat_state = get_num_fields(entity_type);
      auto lower_string = [](const std::string& text) {
        std::string out = text;
        for (char& c : out) { c = char(std::tolower(static_cast<unsigned char>(c))); }
        return out;
      };

      std::vector<int> adjoint_param_field_id(m_element.adjoint_num_params(),-1);
      for (int p=0; p<m_element.adjoint_num_params(); p++) {
        const std::string target_name = "dparam_" + lower_string(m_element.adjoint_param_name(p));
        for (int i=0; i<Nmat_state; i++) {
          if (lower_string(element_field_names[i]) == target_name) {
            adjoint_param_field_id[p] = i;
            break;
          }
        }
      }

      for (int e=0; e<Nelems; e++) {
        Real* elem_state = &state[Nstate*e];
        m_element.m_model.get_fields(elem_state,&field_data[Nmat_state*e]);
        for (int p=0; p<m_element.adjoint_num_params(); p++) {
          if (adjoint_param_field_id[p] >= 0) {
            field_data[Nmat_state*e + adjoint_param_field_id[p]] = m_element.adjoint_get_param_gradient(elem_state,p);
          }
        }
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

  virtual void clear_adjoint_state(void) {
    std::fill(u_adjoint.begin(),u_adjoint.end(),0.0);
    std::fill(v_adjoint.begin(),v_adjoint.end(),0.0);
    std::fill(kick_force_seed.begin(),kick_force_seed.end(),0.0);
    adjoint_mass_damping_gradient = 0.0;
    adjoint_contact_stiffness_gradient = 0.0;
    std::fill(adjoint_param_gradients.begin(),adjoint_param_gradients.end(),0.0);
  }

  // ===================================================================== //

  virtual void add_nodal_velocity_adjoint_seed(const int* node_ids, const double* seed_xy, int num_nodes) {
    if (node_ids == nullptr) {
      std::cout << "add_nodal_velocity_adjoint_seed received null node_ids pointer" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (seed_xy == nullptr) {
      std::cout << "add_nodal_velocity_adjoint_seed received null seed pointer" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (num_nodes < 0) {
      std::cout << "add_nodal_velocity_adjoint_seed received negative node count: " << num_nodes << std::endl;
      exit(EXIT_FAILURE);
    }
    for (int i=0; i<num_nodes; i++) {
      const int node_id = node_ids[i];
      if ((node_id < 0) || (node_id >= Nnodes)) {
        std::cout << "add_nodal_velocity_adjoint_seed received invalid node id: " << node_id << std::endl;
        exit(EXIT_FAILURE);
      }
      const int dof_x = Ndofs_per_node*node_id + 0;
      if ((dof_x >= 0) && (dof_x < Ndofs)) {
        v_adjoint[dof_x] += seed_xy[2*i + 0];
      }
      if (Ndofs_per_node > 1) {
        const int dof_y = Ndofs_per_node*node_id + 1;
        if ((dof_y >= 0) && (dof_y < Ndofs)) {
          v_adjoint[dof_y] += seed_xy[2*i + 1];
        }
      }
    }
  }

  // ===================================================================== //

  virtual void add_nodal_displacement_adjoint_seed(const int* node_ids, const double* seed_xy, int num_nodes) {
    if (node_ids == nullptr) {
      std::cout << "add_nodal_displacement_adjoint_seed received null node_ids pointer" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (seed_xy == nullptr) {
      std::cout << "add_nodal_displacement_adjoint_seed received null seed pointer" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (num_nodes < 0) {
      std::cout << "add_nodal_displacement_adjoint_seed received negative node count: " << num_nodes << std::endl;
      exit(EXIT_FAILURE);
    }
    for (int i=0; i<num_nodes; i++) {
      const int node_id = node_ids[i];
      if ((node_id < 0) || (node_id >= Nnodes)) {
        std::cout << "add_nodal_displacement_adjoint_seed received invalid node id: " << node_id << std::endl;
        exit(EXIT_FAILURE);
      }
      const int dof_x = Ndofs_per_node*node_id + 0;
      if ((dof_x >= 0) && (dof_x < Ndofs)) {
        u_adjoint[dof_x] += seed_xy[2*i + 0];
      }
      if (Ndofs_per_node > 1) {
        const int dof_y = Ndofs_per_node*node_id + 1;
        if ((dof_y >= 0) && (dof_y < Ndofs)) {
          u_adjoint[dof_y] += seed_xy[2*i + 1];
        }
      }
    }
  }

  // ===================================================================== //
private:
  // ===================================================================== //

  void initialize_adjoint_debug_dump_stream() {
    if (m_adjoint_debug_dump_stream.is_open()) {
      m_adjoint_debug_dump_stream.close();
    }
    if (!m_adjoint_debug_dump_enable) { return; }
    if (m_element_stiffness_param_id < 0) {
      std::cout << "adjoint_debug_dump_enable requested, but material does not expose parameter "
                << "'stiffness_scaling_factor'." << std::endl;
      return;
    }

    const std::string filename = "adjoint_stiffness_qp_debug.csv";
    bool has_existing_data = false;
    {
      std::ifstream in(filename);
      has_existing_data = in.good() && (in.peek() != std::ifstream::traits_type::eof());
    }

    m_adjoint_debug_dump_stream.open(filename,std::ios::out | std::ios::app);
    if (!m_adjoint_debug_dump_stream.is_open()) {
      std::cout << "Failed to open adjoint debug dump file: " << filename << std::endl;
      exit(EXIT_FAILURE);
    }
    if (!has_existing_data) {
      m_adjoint_debug_dump_stream
        << "run_id,time_step,time,phase,stage,element_id,qp_id,"
        << "grad_before,grad_after,grad_delta,abs_grad_delta\n";
    }
    m_adjoint_debug_dump_stream.flush();
  }

  void log_adjoint_stiffness_qp_delta(PassPhase phase, const char* stage,
                                      int element_id, int qp_id,
                                      Real grad_before, Real grad_after) {
    if (!m_adjoint_debug_dump_enable) { return; }
    if (m_element_stiffness_param_id < 0) { return; }
    if (!m_adjoint_debug_dump_stream.is_open()) { return; }
    if (m_adjoint_debug_dump_rows >= m_adjoint_debug_dump_max_rows) { return; }
    if ((m_adjoint_debug_dump_stride > 1) && (element_id % m_adjoint_debug_dump_stride != 0)) { return; }

    const Real grad_delta = grad_after - grad_before;
    const Real abs_grad_delta = std::fabs(grad_delta);
    if (abs_grad_delta < m_adjoint_debug_dump_threshold) { return; }

    m_adjoint_debug_dump_stream
      << m_adjoint_debug_dump_run_id << ","
      << m_time_step << ","
      << std::setprecision(17) << m_time << ","
      << pass_phase_name(phase) << ","
      << stage << ","
      << element_id << ","
      << qp_id << ","
      << std::setprecision(17) << grad_before << ","
      << std::setprecision(17) << grad_after << ","
      << std::setprecision(17) << grad_delta << ","
      << std::setprecision(17) << abs_grad_delta << "\n";
    m_adjoint_debug_dump_rows++;

    if (m_adjoint_debug_dump_rows == m_adjoint_debug_dump_max_rows) {
      std::cout << "adjoint debug dump reached max rows (" << m_adjoint_debug_dump_max_rows
                << ") for run_id=" << m_adjoint_debug_dump_run_id << std::endl;
    }
  }

  LongInteger dual_mantissa_abs(const Real&) const { return 0; }

  template<typename T>
  LongInteger dual_mantissa_abs(const T& value) const {
    return std::llabs(LongInteger(value.mantissa));
  }

  void maybe_warn_dual_velocity_overflow(const char* context) {
    if (!m_dual_overflow_warn_enable) { return; }
    if (m_dual_overflow_warn_limit == 0) { return; }
    if (m_dual_overflow_warn_count >= m_dual_overflow_warn_limit) { return; }

    const LongInteger max_integer = LongInteger(std::numeric_limits<Integer>::max());
    const LongInteger threshold = LongInteger(m_dual_overflow_warn_fraction*Real(max_integer));

    LongInteger max_abs_mantissa = 0;
    int max_abs_dof = -1;
    for (int i=0; i<Ndofs; i++) {
      LongInteger abs_value = dual_mantissa_abs(v[i].second);
      if (abs_value > max_abs_mantissa) {
        max_abs_mantissa = abs_value;
        max_abs_dof = i;
      }
    }

    if (max_abs_mantissa >= threshold) {
      m_dual_overflow_warn_count++;
      std::cout
        << "WARNING: dual velocity mantissa near overflow in context=" << context
        << " at time_step=" << m_time_step
        << ", time=" << m_time
        << ", dof=" << max_abs_dof
        << ", |mantissa|max=" << max_abs_mantissa
        << ", threshold=" << threshold
        << ", int32_max=" << max_integer
        << std::endl;
      if (m_dual_overflow_warn_count == m_dual_overflow_warn_limit) {
        std::cout << "Further dual-overflow warnings suppressed for this run." << std::endl;
      }
    }
  }

  void register_adjoint_param_name(const std::string& name) {
    if (name.empty()) { return; }
    if (std::find(adjoint_param_names.begin(),adjoint_param_names.end(),name) == adjoint_param_names.end()) {
      adjoint_param_names.push_back(name);
      adjoint_param_gradients.push_back(0.0);
      global_field_names.push_back("dL_dparam_" + name);
    }
  }

  int adjoint_param_index(const std::string& name) const {
    for (int i=0; i<int(adjoint_param_names.size()); i++) {
      if (adjoint_param_names[i] == name) { return i; }
    }
    return -1;
  }

  void register_adjoint_params_from_materials() {
    register_adjoint_param_name("mass_damping_factor");
    register_adjoint_param_name("contact_stiffness");
    for (int i=0; i<m_element.adjoint_num_params(); i++) {
      register_adjoint_param_name(m_element.adjoint_param_name(i));
    }
    for (int i=0; i<m_truss.adjoint_num_params(); i++) {
      register_adjoint_param_name(m_truss.adjoint_param_name(i));
    }
  }

  void update_global_param_gradients() {
    std::fill(adjoint_param_gradients.begin(),adjoint_param_gradients.end(),0.0);

    const int Nstate_vars_per_elem = m_element.num_state_vars();
    for (int e=0; e<Nelems; e++) {
      Real* elem_state = &state[Nstate_vars_per_elem*e];
      for (int p=0; p<m_element.adjoint_num_params(); p++) {
        const std::string pname = m_element.adjoint_param_name(p);
        int idx = adjoint_param_index(pname);
        if (idx >= 0) {
          adjoint_param_gradients[idx] += m_element.adjoint_get_param_gradient(elem_state,p);
        }
      }
    }

    const int Nstate_vars_per_truss = m_truss.num_state_vars();
    for (int e=0; e<Ntruss; e++) {
      Real* truss_state_e = &truss_state[Nstate_vars_per_truss*e];
      for (int p=0; p<m_truss.adjoint_num_params(); p++) {
        const std::string pname = m_truss.adjoint_param_name(p);
        int idx = adjoint_param_index(pname);
        if (idx >= 0) {
          adjoint_param_gradients[idx] += m_truss.adjoint_get_param_gradient(truss_state_e,p);
        }
      }
    }

    int alpha_idx = adjoint_param_index("mass_damping_factor");
    if (alpha_idx >= 0) {
      adjoint_param_gradients[alpha_idx] += adjoint_mass_damping_gradient;
    }
    int wall_idx = adjoint_param_index("contact_stiffness");
    if (wall_idx >= 0) {
      adjoint_param_gradients[wall_idx] += adjoint_contact_stiffness_gradient;
    }
  }

  void enforce_backward_adjoint_guardrails() const {
    if (!m_contact_interactions.empty()) {
      std::cout << "BackwardAdjoint currently supports only internal-force adjoint terms; "
                << "segment-contact interaction adjoint terms are not included in Slice D." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void apply_second_kick_damping_pullback(Real dt_step) {
    if ((dt_step <= 0.0) || (std::fabs(m_alpha) == 0.0)) { return; }
    const Real half_dt = 0.5*dt_step;
    for (int i=0; i<Ndofs; i++) {
      if (fixity[i]) { continue; }
      if (has_time_varying_bc(i)) { continue; }
      const Real phi2 = 1.0 + half_dt*alpha[i];
      if (phi2 == 0.0) {
        std::cout << "BackwardAdjoint damping pullback encountered zero phi2 at dof " << i << std::endl;
        exit(EXIT_FAILURE);
      }
      const Real bar_v_np1 = v_adjoint[i];
      const Real bar_v_half = bar_v_np1/phi2;
      adjoint_mass_damping_gradient += -half_dt*bar_v_half*Real(v[i].first);
      v_adjoint[i] = bar_v_half;
    }
  }

  void apply_first_kick_damping_pullback(Real dt_step) {
    if ((dt_step <= 0.0) || (std::fabs(m_alpha) == 0.0)) { return; }
    const Real half_dt = 0.5*dt_step;
    for (int i=0; i<Ndofs; i++) {
      if (fixity[i]) { continue; }
      if (has_time_varying_bc(i)) { continue; }
      const Real phi1 = 1.0 - half_dt*alpha[i];
      const Real bar_v_half = v_adjoint[i];
      adjoint_mass_damping_gradient += -half_dt*bar_v_half*Real(v[i].first);
      v_adjoint[i] = phi1*bar_v_half;
    }
  }

  void build_kick_force_seed(Real dt_step, const std::vector<Real>& mass_values) {
    std::fill(kick_force_seed.begin(),kick_force_seed.end(),0.0);
    if (dt_step <= 0.0) { return; }
    for (int i=0; i<Ndofs; i++) {
      if (fixity[i]) { continue; }
      if (has_time_varying_bc(i)) { continue; }
      if (mass_values[i] <= 0.0) { continue; }
      kick_force_seed[i] = 0.5*dt_step*v_adjoint[i]/mass_values[i];
    }
  }

  void compute_mass_vector(std::vector<Real>& mass_values) {
    mass_values.assign(Ndofs,0.0);

    // Solid element mass contributions (same midpoint rule as Element::update).
    for (int e=0; e<Nelems; e++) {
      Real xe[8];
      for (int j=0; j<Nnodes_per_elem; j++) {
        const int jnode_id = connect[Nnodes_per_elem*e+j];
        for (int i=0; i<Ndofs_per_node; i++) {
          xe[Ndofs_per_node*j+i] = x[Ndofs_per_node*jnode_id+i];
        }
      }

      Real dxi[4]  = { -0.25, +0.25, +0.25, -0.25 };
      Real deta[4] = { -0.25, -0.25, +0.25, +0.25 };
      Real J[2][2] = { {0.0,0.0}, {0.0,0.0} };
      for (int i=0; i<4; i++) {
        J[0][0] += xe[2*i+0] * dxi[i]; J[0][1] += xe[2*i+0] * deta[i];
        J[1][0] += xe[2*i+1] * dxi[i]; J[1][1] += xe[2*i+1] * deta[i];
      }
      Real detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
      Real mass = detJ*m_element.m_model.density();
      for (int j=0; j<4; j++) {
        const int jnode_id = connect[Nnodes_per_elem*e+j];
        mass_values[Ndofs_per_node*jnode_id+0] += mass;
        mass_values[Ndofs_per_node*jnode_id+1] += mass;
      }
    }

    // Truss element mass contributions (same midpoint rule as Truss::update).
    for (int e=0; e<Ntruss; e++) {
      const int n0 = truss_connect[2*e+0];
      const int n1 = truss_connect[2*e+1];
      Real dxi[2] = { -0.5, +0.5 };
      Real J[2] = { 0.0, 0.0 };
      J[0] += x[Ndofs_per_node*n0+0] * dxi[0];
      J[1] += x[Ndofs_per_node*n0+1] * dxi[0];
      J[0] += x[Ndofs_per_node*n1+0] * dxi[1];
      J[1] += x[Ndofs_per_node*n1+1] * dxi[1];
      Real normJ = std::sqrt(J[0]*J[0] + J[1]*J[1]);
      Real mass = normJ*m_truss.m_model.mass_per_unit_length();
      mass_values[Ndofs_per_node*n0+0] += mass;
      mass_values[Ndofs_per_node*n0+1] += mass;
      mass_values[Ndofs_per_node*n1+0] += mass;
      mass_values[Ndofs_per_node*n1+1] += mass;
    }

    for (int i=0; i<Npoints; i++) {
      mass_values[2*point_ids[i]+0] += point_mass[i];
      mass_values[2*point_ids[i]+1] += point_mass[i];
    }
  }

  void apply_global_kick_adjoint(Real dt_step, bool first_kick) {
    if (!first_kick) {
      build_kick_force_seed(dt_step,m);
    }

    // Early out when there is no kick seed.
    Real max_seed = 0.0;
    for (int i=0; i<Ndofs; i++) max_seed = std::max(max_seed,std::fabs(kick_force_seed[i]));
    if (max_seed == 0.0) { return; }

    // 2D element contribution.
    const int Nstate_vars_per_elem = m_element.num_state_vars();
    const int Nmat_state = m_element.m_model.num_state_vars();
    for (int e=0; e<Nelems; e++) {
      Real xe[8];
      Real ue[8];
      int node_ids[4];
      for (int j=0; j<Nnodes_per_elem; j++) {
        const int jnode_id = connect[Nnodes_per_elem*e+j];
        node_ids[j] = jnode_id;
        xe[2*j+0] = x[Ndofs_per_node*jnode_id+0];
        xe[2*j+1] = x[Ndofs_per_node*jnode_id+1];
        ue[2*j+0] = u[Ndofs_per_node*jnode_id+0].first;
        ue[2*j+1] = u[Ndofs_per_node*jnode_id+1].first;
      }

      const Real sqrt_third = 1.0/std::sqrt(3.0);
      Real xi[4] = { -1.0, +1.0, +1.0, -1.0 };
      Real eta[4] = { -1.0, -1.0, +1.0, +1.0 };

      for (int q=0; q<4; q++) {
        Real xiq = xi[q]*sqrt_third;
        Real etaq = eta[q]*sqrt_third;

        Real dxi[4], deta[4];
        for (int i=0; i<4; i++) {
          dxi[i]  = 0.25*xi[i]*(1.0+etaq*eta[i]);
          deta[i] = 0.25*eta[i]*(1.0+xiq*xi[i]);
        }

        Real J0[2][2] = { {0.0,0.0}, {0.0,0.0} };
        Real J[2][2]  = { {0.0,0.0}, {0.0,0.0} };
        for (int i=0; i<4; i++) {
          J0[0][0] += xe[2*i+0] * dxi[i]; J0[0][1] += xe[2*i+0] * deta[i];
          J0[1][0] += xe[2*i+1] * dxi[i]; J0[1][1] += xe[2*i+1] * deta[i];
          Real xtx = xe[2*i+0] + ue[2*i+0];
          Real xty = xe[2*i+1] + ue[2*i+1];
          J[0][0] += xtx * dxi[i]; J[0][1] += xtx * deta[i];
          J[1][0] += xty * dxi[i]; J[1][1] += xty * deta[i];
        }

        Real inv_detJ0 = 1.0/(J0[0][0]*J0[1][1] - J0[0][1]*J0[1][0]);
        Real invJ0[2][2];
        invJ0[0][0] = +J0[1][1]*inv_detJ0;
        invJ0[0][1] = -J0[0][1]*inv_detJ0;
        invJ0[1][0] = -J0[1][0]*inv_detJ0;
        invJ0[1][1] = +J0[0][0]*inv_detJ0;

        Real* model_state = &state[Nstate_vars_per_elem*e + (q+1)*Nmat_state];
        Real sxx = model_state[0];
        Real syy = model_state[1];
        Real sxy = model_state[5];

        Real bar_P00 = 0.0, bar_P01 = 0.0, bar_P10 = 0.0, bar_P11 = 0.0;
        for (int i=0; i<4; i++) {
          Real bar_fe_x = -kick_force_seed[Ndofs_per_node*node_ids[i]+0];
          Real bar_fe_y = -kick_force_seed[Ndofs_per_node*node_ids[i]+1];
          bar_P00 += bar_fe_x*dxi[i];
          bar_P01 += bar_fe_x*deta[i];
          bar_P10 += bar_fe_y*dxi[i];
          bar_P11 += bar_fe_y*deta[i];
        }

        Real cof00 = +J[1][1];
        Real cof10 = -J[0][1];
        Real cof01 = -J[1][0];
        Real cof11 = +J[0][0];

        Real bar_sxx = bar_P00*cof00 + bar_P01*cof01;
        Real bar_syy = bar_P10*cof10 + bar_P11*cof11;
        Real bar_sxy = bar_P00*cof10 + bar_P01*cof11 + bar_P10*cof00 + bar_P11*cof01;

        if (!first_kick) {
          Real grad_before_direct = 0.0;
          if (m_adjoint_debug_dump_enable && (m_element_stiffness_param_id >= 0)) {
            grad_before_direct =
              material_adjoint_get_param_gradient(m_element.m_model,model_state,m_element_stiffness_param_id);
          }
          material_adjoint_add_history_seed_from_stress(m_element.m_model,model_state,bar_sxx,bar_syy,bar_sxy);
          material_adjoint_add_direct_param_seed_from_stress(m_element.m_model,model_state,bar_sxx,bar_syy,bar_sxy);
          if (m_adjoint_debug_dump_enable && (m_element_stiffness_param_id >= 0)) {
            Real grad_after_direct =
              material_adjoint_get_param_gradient(m_element.m_model,model_state,m_element_stiffness_param_id);
            log_adjoint_stiffness_qp_delta(
              PassPhase::BackwardAdjoint,
              "second_kick_direct",
              e,
              q,
              grad_before_direct,
              grad_after_direct
            );
          }
        }

        Real bar_cof00 = bar_P00*sxx + bar_P10*sxy;
        Real bar_cof01 = bar_P01*sxx + bar_P11*sxy;
        Real bar_cof10 = bar_P00*sxy + bar_P10*syy;
        Real bar_cof11 = bar_P01*sxy + bar_P11*syy;

        Real bar_exx = 0.0, bar_eyy = 0.0, bar_gxy = 0.0;
        material_adjoint_pullback_stress_to_strain(m_element.m_model,model_state,
                                                   bar_sxx,bar_syy,bar_sxy,
                                                   bar_exx,bar_eyy,bar_gxy);

        Real bar_F00 = bar_exx;
        Real bar_F11 = bar_eyy;
        Real bar_F01 = bar_gxy;
        Real bar_F10 = bar_gxy;

        Real bar_J00 = bar_F00*invJ0[0][0] + bar_F01*invJ0[0][1];
        Real bar_J01 = bar_F00*invJ0[1][0] + bar_F01*invJ0[1][1];
        Real bar_J10 = bar_F10*invJ0[0][0] + bar_F11*invJ0[0][1];
        Real bar_J11 = bar_F10*invJ0[1][0] + bar_F11*invJ0[1][1];

        // cof(J) pullback
        bar_J11 += bar_cof00;
        bar_J01 += -bar_cof10;
        bar_J10 += -bar_cof01;
        bar_J00 += bar_cof11;

        for (int i=0; i<4; i++) {
          u_adjoint[Ndofs_per_node*node_ids[i]+0] += bar_J00*dxi[i] + bar_J01*deta[i];
          u_adjoint[Ndofs_per_node*node_ids[i]+1] += bar_J10*dxi[i] + bar_J11*deta[i];
        }
      }
    }

    // Truss contribution.
    const int Nstate_vars_per_truss = m_truss.num_state_vars();
    for (int e=0; e<Ntruss; e++) {
      int n0 = truss_connect[2*e+0];
      int n1 = truss_connect[2*e+1];
      Real xe[4] = { x[2*n0+0], x[2*n0+1], x[2*n1+0], x[2*n1+1] };
      Real ue[4] = { u[2*n0+0].first, u[2*n0+1].first, u[2*n1+0].first, u[2*n1+1].first };
      Real* state_e = &truss_state[Nstate_vars_per_truss*e];

      const Real dxi[2] = { -0.5, +0.5 };
      const Real wq = 2.0;
      const Real c0 = dxi[0]*wq;
      const Real c1 = dxi[1]*wq;

      Real J0x = xe[0]*dxi[0] + xe[2]*dxi[1];
      Real J0y = xe[1]*dxi[0] + xe[3]*dxi[1];
      Real Jx = (xe[0]+ue[0])*dxi[0] + (xe[2]+ue[2])*dxi[1];
      Real Jy = (xe[1]+ue[1])*dxi[0] + (xe[3]+ue[3])*dxi[1];
      Real normJ0 = std::sqrt(J0x*J0x + J0y*J0y);
      Real normJ = std::sqrt(Jx*Jx + Jy*Jy);
      if ((normJ0 <= 0.0) || (normJ <= 0.0)) { continue; }

      Real tx = Jx/normJ;
      Real ty = Jy/normJ;
      Real sigma = state_e[0];

      Real bar_fe0x = -kick_force_seed[2*n0+0];
      Real bar_fe0y = -kick_force_seed[2*n0+1];
      Real bar_fe1x = -kick_force_seed[2*n1+0];
      Real bar_fe1y = -kick_force_seed[2*n1+1];

      Real bar_sigma = c0*(bar_fe0x*tx + bar_fe0y*ty) + c1*(bar_fe1x*tx + bar_fe1y*ty);
      if (!first_kick) {
        material_adjoint_add_history_seed_from_stress(m_truss.m_model,state_e,bar_sigma);
        material_adjoint_add_direct_param_seed_from_stress(m_truss.m_model,state_e,bar_sigma);
      }

      Real bar_strain = 0.0;
      material_adjoint_pullback_stress_to_strain(m_truss.m_model,state_e,bar_sigma,bar_strain);

      Real bar_lambda = bar_strain;
      Real bar_normJ = bar_lambda/normJ0;
      Real bar_Jx = bar_normJ*tx;
      Real bar_Jy = bar_normJ*ty;

      Real bar_tx = sigma*(c0*bar_fe0x + c1*bar_fe1x);
      Real bar_ty = sigma*(c0*bar_fe0y + c1*bar_fe1y);
      Real bdot = bar_tx*tx + bar_ty*ty;
      bar_Jx += (bar_tx - bdot*tx)/normJ;
      bar_Jy += (bar_ty - bdot*ty)/normJ;

      u_adjoint[2*n0+0] += dxi[0]*bar_Jx;
      u_adjoint[2*n0+1] += dxi[0]*bar_Jy;
      u_adjoint[2*n1+0] += dxi[1]*bar_Jx;
      u_adjoint[2*n1+1] += dxi[1]*bar_Jy;
    }

    // Rigid-wall contact contribution (f_y += k * max(0,-y)).
    if ((Ndofs_per_node > 1) && (std::fabs(m_contact_stiffness) > 0.0)) {
      for (int i=0; i<Nnodes; i++) {
        const int dof_y = Ndofs_per_node*i+1;
        const Real y = x[dof_y] + Real(u[dof_y].first);
        if (y < 0.0) {
          const Real penetration = -y;
          const Real bar_fy = kick_force_seed[dof_y];
          u_adjoint[dof_y] += -m_contact_stiffness*bar_fy;
          adjoint_contact_stiffness_gradient += penetration*bar_fy;
        }
      }
    }
  }

  // Check if a specified DoF is controlled by a time-varying boundary condition
  bool has_time_varying_bc(int dof) const { return (dof >= 0 && dof < int(m_has_time_bc.size())) ? m_has_time_bc[dof] : false; }

  // ===================================================================== //

  // Apply all time-varying displacement boundary conditions at the specified time
  void apply_displacement_bcs(Real target_time, Real dt) {
    if (m_displacement_bcs.empty()) { return; }
    
    // later if we want to increase the number of elements/nodes we should try to save num_nodes to size_t so it will be size_t i and j for loops or even dofs
    for (int i = 0; i < m_displacement_bcs.size(); i++) {
      DisplacementBoundaryCondition& bc = m_displacement_bcs[i];
      if (bc.function == nullptr) { continue; }
      int num_nodes = bc.nodes.node_ids.size();
      if (bc.last_values.size() != num_nodes) {
	bc.last_values.resize(num_nodes,0.0);
      }
      if (bc.prescribed_velocities.size() != num_nodes) {
	bc.prescribed_velocities.resize(num_nodes,0.0);
      }
      // later if we want to increase the number of elements/nodes we should try to save num_nodes to size_t so it will be size_t i and j for loops or even dofs
      for (int j=0; j<num_nodes; ++j) {
	int node_id = bc.nodes.node_ids[j];
	if ((node_id < 0) || (node_id >= Nnodes)) { continue; }

	int dof = Ndofs_per_node*node_id + bc.component;
	if ((dof < 0) || (dof >= Ndofs)) { continue; }

	Real px = x[Ndofs_per_node*node_id + 0];
	Real py = 0.0;
	if (Ndofs_per_node > 1) py = x[Ndofs_per_node*node_id + 1];

	Real previous_diplacement_value = bc.last_values[j];
	Real new_diplacement_value = bc.function(target_time, px, py);

	Real velocity = 0.0;
	if (dt != 0.0) velocity = (new_diplacement_value - previous_diplacement_value)/dt;
	bc.last_values[j] = new_diplacement_value;
	bc.prescribed_velocities[j] = velocity;
	u[dof] = Dual<FixedU>(new_diplacement_value,0.0);
	v[dof] = Dual<FixedV>(velocity,0.0);
	m_has_time_bc[dof] = true;
	fixity[dof] = true;
      }
    }
  } // apply_displacement_bcs()

  // ===================================================================== //

  void enforce_prescribed_velocities() {
    if (m_displacement_bcs.empty()) { return; }

    for (int i = 0; i < m_displacement_bcs.size(); i++) {
      DisplacementBoundaryCondition& bc = m_displacement_bcs[i];
      int num_nodes = bc.nodes.node_ids.size();
      for (int j=0; j<num_nodes; j++) {
        int node_id = bc.nodes.node_ids[j];
        if ((node_id < 0) || (node_id >= Nnodes)) continue; 
        int dof = Ndofs_per_node*node_id + bc.component;
        if ((dof < 0) || (dof >= Ndofs)) continue; 
        if (!has_time_varying_bc(dof)) continue; 
        if (j < bc.prescribed_velocities.size()) {
          v[dof].first = bc.prescribed_velocities[j];
          v[dof].second = 0.0;
        }
      }
    }
  } // enforce_prescribed_velocities()

  // ===================================================================== //

  // Procedure to update the accelerations
  void update_accelerations(Real dt, PassPhase phase, Real dt_step) {
    const bool first_kick_adjoint = (phase == PassPhase::BackwardAdjoint) && (dt < 0.0);
    if (first_kick_adjoint) {
      std::vector<Real> kick_mass;
      compute_mass_vector(kick_mass);
      build_kick_force_seed(dt_step,kick_mass);
    } else {
      std::fill(kick_force_seed.begin(),kick_force_seed.end(),0.0);
    }

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
    const int Nmat_state_elem = m_element.m_model.num_state_vars();

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
      m_element.adjoint_clear_step_seed(&state[Nstate_vars_per_elem*e]);
      if (phase == PassPhase::BackwardAdjoint) {
        m_element.adjoint_objective_seed(&state[Nstate_vars_per_elem*e]);
        if (first_kick_adjoint) {
          const int Nmat_state = m_element.m_model.num_state_vars();
          const Real sqrt_third = 1.0/std::sqrt(3.0);
          Real xi[4] = { -1.0, +1.0, +1.0, -1.0 };
          Real eta[4] = { -1.0, -1.0, +1.0, +1.0 };
          for (int q=0; q<4; q++) {
            Real xiq = xi[q]*sqrt_third;
            Real etaq = eta[q]*sqrt_third;

            Real dxi[4], deta[4];
            for (int i=0; i<4; i++) {
              dxi[i]  = 0.25*xi[i]*(1.0+etaq*eta[i]);
              deta[i] = 0.25*eta[i]*(1.0+xiq*xi[i]);
            }

            Real J0[2][2] = { {0.0,0.0}, {0.0,0.0} };
            Real J[2][2]  = { {0.0,0.0}, {0.0,0.0} };
            int node_ids[4];
            for (int i=0; i<4; i++) {
              node_ids[i] = connect[Nnodes_per_elem*e+i];
              J0[0][0] += xe[2*i+0] * dxi[i]; J0[0][1] += xe[2*i+0] * deta[i];
              J0[1][0] += xe[2*i+1] * dxi[i]; J0[1][1] += xe[2*i+1] * deta[i];
              Real xtx = xe[2*i+0] + ue[2*i+0];
              Real xty = xe[2*i+1] + ue[2*i+1];
              J[0][0] += xtx * dxi[i]; J[0][1] += xtx * deta[i];
              J[1][0] += xty * dxi[i]; J[1][1] += xty * deta[i];
            }

            Real bar_P00 = 0.0, bar_P01 = 0.0, bar_P10 = 0.0, bar_P11 = 0.0;
            for (int i=0; i<4; i++) {
              Real bar_fe_x = -kick_force_seed[Ndofs_per_node*node_ids[i]+0];
              Real bar_fe_y = -kick_force_seed[Ndofs_per_node*node_ids[i]+1];
              bar_P00 += bar_fe_x*dxi[i];
              bar_P01 += bar_fe_x*deta[i];
              bar_P10 += bar_fe_y*dxi[i];
              bar_P11 += bar_fe_y*deta[i];
            }

            Real cof00 = +J[1][1];
            Real cof10 = -J[0][1];
            Real cof01 = -J[1][0];
            Real cof11 = +J[0][0];

            Real bar_sxx = bar_P00*cof00 + bar_P01*cof01;
            Real bar_syy = bar_P10*cof10 + bar_P11*cof11;
            Real bar_sxy = bar_P00*cof10 + bar_P01*cof11 + bar_P10*cof00 + bar_P11*cof01;

            Real* model_state = &state[Nstate_vars_per_elem*e + (q+1)*Nmat_state];
            material_adjoint_add_stress_seed(m_element.m_model,model_state,bar_sxx,bar_syy,bar_sxy);
          }
        }
      }
      Real stiffness_grad_before_update[4] = { 0.0, 0.0, 0.0, 0.0 };
      if ((phase == PassPhase::BackwardAdjoint) && m_adjoint_debug_dump_enable && (m_element_stiffness_param_id >= 0)) {
        for (int q=0; q<4; q++) {
          Real* model_state = &state[Nstate_vars_per_elem*e + (q+1)*Nmat_state_elem];
          stiffness_grad_before_update[q] =
            material_adjoint_get_param_gradient(m_element.m_model,model_state,m_element_stiffness_param_id);
        }
      }

      m_element.update(xe,ue,me,fe,Ee,&state[Nstate_vars_per_elem*e],dt_step,phase);

      if ((phase == PassPhase::BackwardAdjoint) && m_adjoint_debug_dump_enable && (m_element_stiffness_param_id >= 0)) {
        for (int q=0; q<4; q++) {
          Real* model_state = &state[Nstate_vars_per_elem*e + (q+1)*Nmat_state_elem];
          Real stiffness_grad_after_update =
            material_adjoint_get_param_gradient(m_element.m_model,model_state,m_element_stiffness_param_id);
          log_adjoint_stiffness_qp_delta(
            phase,
            "material_update",
            e,
            q,
            stiffness_grad_before_update[q],
            stiffness_grad_after_update
          );
        }
      }

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
      m_truss.adjoint_clear_step_seed(&truss_state[Nstate_vars_per_truss*e]);
      if (phase == PassPhase::BackwardAdjoint) {
        m_truss.adjoint_objective_seed(&truss_state[Nstate_vars_per_truss*e]);
        if (first_kick_adjoint) {
          const int n0 = truss_connect[2*e+0];
          const int n1 = truss_connect[2*e+1];
          const Real dxi[2] = { -0.5, +0.5 };
          const Real wq = 2.0;
          const Real c0 = dxi[0]*wq;
          const Real c1 = dxi[1]*wq;

          Real Jx = (xe[0]+ue[0])*dxi[0] + (xe[2]+ue[2])*dxi[1];
          Real Jy = (xe[1]+ue[1])*dxi[0] + (xe[3]+ue[3])*dxi[1];
          Real normJ = std::sqrt(Jx*Jx + Jy*Jy);
          if (normJ > 0.0) {
            Real tx = Jx/normJ;
            Real ty = Jy/normJ;
            Real bar_fe0x = -kick_force_seed[2*n0+0];
            Real bar_fe0y = -kick_force_seed[2*n0+1];
            Real bar_fe1x = -kick_force_seed[2*n1+0];
            Real bar_fe1y = -kick_force_seed[2*n1+1];
            Real bar_sigma = c0*(bar_fe0x*tx + bar_fe0y*ty) + c1*(bar_fe1x*tx + bar_fe1y*ty);
            material_adjoint_add_stress_seed(m_truss.m_model,&truss_state[Nstate_vars_per_truss*e],bar_sigma);
          }
        }
      }
      m_truss.update(xe,ue,me,fe,Ee,&truss_state[Nstate_vars_per_truss*e],dt_step,phase);

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

    if (first_kick_adjoint) {
      // First kick adjoint at reconstructed state n:
      // material seeds were injected before local backward updates;
      // this call adds global geometric/kinematic pullback terms.
      apply_global_kick_adjoint(dt_step,true);
    }

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

    maybe_warn_dual_velocity_overflow("update_accelerations");
    
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
