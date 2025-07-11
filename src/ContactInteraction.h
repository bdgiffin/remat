#ifndef CONTACT_INTERACTION_H
#define CONTACT_INTERACTION_H

#include "types.h"
#include "Parameters.h"
#include <vector>
#include <map>
#include <limits>
#include <algorithm> // For std::sort, std::unique, std::find, std::fill
#include <math.h>

struct ContactInteraction {

  // Nodal contact stiffness parameter
  Real m_contact_stiffness = 0.0;
  Real m_search_radius     = std::numeric_limits<Real>::max();

  // The global node IDs of all nodes in the contact interaction
  int Nnodes = 0;
  std::vector<int> nodes;
  std::vector<Real> x_nodes; // current coordinates of all nodes
  std::vector<Real> f_nodes; // accumulated forces of all nodes

  // The segment connectivities of all surfaces in the contact interaction
  int Nsegments = 0;
  std::vector<int> segments;

  // The nodes of all segments
  int Nsegment_nodes = 0;
  std::vector<int> segment_nodes; // sorted list of unique segment node global IDs
  std::vector<int> node_segments; // list of (up to 2) segments associated with each unique node
  std::vector<Real> x_segment_nodes; // current coordinates of all segment_nodes
  std::vector<Real> f_segment_nodes; // accumulated forces of all segment_nodes

  // Empty constructor
  ContactInteraction(void) { }

  // Parameterized constructor
  ContactInteraction(Parameters& params) {
    if (params.count("contact_stiffness") > 0) m_contact_stiffness = params["contact_stiffness"];
    if (params.count("search_radius") > 0)     m_search_radius     = params["search_radius"];
  }

  // Initialize contact interaction with list of nodes and segment connectivity
  void initialize(int* node_ids, int* segment_connectivity, int new_Nnodes, int new_Nsegments, Parameters& params) {
    // Set contact stiffness and search radius parameters, if not already set
    if (params.count("contact_stiffness") > 0) m_contact_stiffness = params["contact_stiffness"];
    if (params.count("search_radius") > 0)     m_search_radius     = params["search_radius"];

    // Set totals
    Nnodes    = new_Nnodes;
    Nsegments = new_Nsegments;
    
    // Initialize arrays
    nodes.resize(Nnodes);
    segments.resize(2*Nsegments);
    for (int i=0; i<Nnodes; i++) nodes[i] = node_ids[i];
    for (int i=0; i<Nsegments; i++) {
      segments[2*i+0] = segment_connectivity[2*i+0];
      segments[2*i+1] = segment_connectivity[2*i+1];
    }

    // Find unique segment node IDs
    segment_nodes = segments;
    std::sort(segment_nodes.begin(), segment_nodes.end());
    auto newEnd = std::unique(segment_nodes.begin(), segment_nodes.end());
    segment_nodes.erase(newEnd, segment_nodes.end());
    Nsegment_nodes = segment_nodes.size();
    
    // Determine the segments associated with each node
    node_segments = std::vector<int>(2*Nsegment_nodes,-1);
    for (int i=0; i<Nsegments; i++) {
      for (int j=0; j<2; j++) {
	int jnode = segments[2*i+j];
	auto it = std::find(segment_nodes.begin(), segment_nodes.end(), jnode);
	int index = std::distance(segment_nodes.begin(), it);
	int k = 0;
	if (node_segments[2*index+k] != -1) k++;
	node_segments[2*index+k] = i;

	// define "segments" array in terms of the locally sorted segment_node indices
	segments[2*i+j] = index;
      }
    }

    // Allocate coordinates and forces for all nodes and segment_nodes
    x_nodes.resize(2*Nnodes);
    f_nodes.resize(2*Nnodes);
    x_segment_nodes.resize(2*Nsegment_nodes);
    f_segment_nodes.resize(2*Nsegment_nodes);
    
  } // initialize()

  // Update nodal force contributions based on the current deformed coordinates of all nodes
  void update_contact_forces(Real* x, Real* f, Real dt) {
    const int Ndofs_per_node = 2;
    
    // Zero-initialize contact forces
    std::fill(f_nodes.begin(),         f_nodes.end(),         0.0);
    std::fill(f_segment_nodes.begin(), f_segment_nodes.end(), 0.0);
    
    // Copy the updated coordinates for all nodes and segment_nodes
    for (int i=0; i<Nnodes; i++) {
      x_nodes[Ndofs_per_node*i+0] = x[Ndofs_per_node*nodes[i]+0];
      x_nodes[Ndofs_per_node*i+1] = x[Ndofs_per_node*nodes[i]+1];
    }
    for (int i=0; i<Nsegment_nodes; i++) {
      x_segment_nodes[Ndofs_per_node*i+0] = x[Ndofs_per_node*segment_nodes[i]+0];
      x_segment_nodes[Ndofs_per_node*i+1] = x[Ndofs_per_node*segment_nodes[i]+1];
    }

    // Resolve contact interactions for all nodes
    for (int j=0; j<nodes.size(); j++) {
      // Get the coordinates of the current node
      Real x0 = x_nodes[Ndofs_per_node*j+0];
      Real y0 = x_nodes[Ndofs_per_node*j+1];
      
      // Find the nearest segment_node to each node (by an exhaustive search)
      int nearest;
      Real min_distance = std::numeric_limits<Real>::max();
      for (int i = 0; i < segment_nodes.size(); i++) {
	Real dx = x0 - x_segment_nodes[Ndofs_per_node*i+0];
	Real dy = y0 - x_segment_nodes[Ndofs_per_node*i+1];
	Real distance = std::sqrt(dx*dx + dy*dy);
	if (distance < min_distance) {
	  min_distance = distance;
	  nearest = i;
	}
      }

      // Only apply forces if the nearest node lies within the specified search radius
      if (min_distance < m_search_radius) {

	// Loop over adjoining segments
	for (int s=0; s<2; s++) {
	  // get the segment ID of the current segment
	  int sid = node_segments[2*nearest+s];
	  if (sid != -1) { // check if the segment exists
	    // get the indices of the segment nodes
	    int i1 = segments[2*sid+0]; // index of segment node 1
	    int i2 = segments[2*sid+1]; // index of segment node 2
	    
	    // get the segment node coordinates (shifted by the node coordinate)
	    Real x1 = x_segment_nodes[2*i1+0] - x0;
	    Real y1 = x_segment_nodes[2*i1+0] - y0;
	    Real x2 = x_segment_nodes[2*i2+0] - x0;
	    Real y2 = x_segment_nodes[2*i2+0] - y0;

	    // determine the segment tangent and length
	    Real dx = x2 - x1;
	    Real dy = y2 - y1;
	    Real invL2 = 1.0/(dx*dx + dy*dy);
	    Real tx = +dx*invL2;
	    Real ty = +dy*invL2;

	    // determine the normal vector and gap
	    Real invL = std::sqrt(invL2);
	    Real nx = +dy*invL;
	    Real ny = -dx*invL;
	    Real g  = -(x1*nx + y1*ny);

	    // determine the projected coordinates on the segment
	    Real xi1 = 1.0 + x1*tx + y1*ty;
	    Real xi2 = 1.0 - xi1;

	    // only enforce contact if the node's projected location lies on the segment and the gap is negative
	    if ((xi1 >= 0.0) && (xi1 <= 1.0) && (g < 0.0)) {
	      Real fn = m_contact_stiffness*g; // normal contact force
	      Real fx = fn*nx; // x-component of contact force applied to node j
	      Real fy = fn*ny; // y-component of contact force applied to node j

	      // sum force contributions to node and corresponding segment nodes
	      f_nodes[2*j+0] += fx;
	      f_nodes[2*j+1] += fy;
	      f_segment_nodes[2*i1+0] -= xi1*fx;
	      f_segment_nodes[2*i1+1] -= xi1*fy;
	      f_segment_nodes[2*i2+0] -= xi2*fx;
	      f_segment_nodes[2*i2+1] -= xi2*fy;
	    }
	    
	  } // if (sid != -1)
	} // for (int s=0; s<2; s++)

      } // if (min_distance < m_search_radius)
      
    } // for (int j=0; j<nodes.size(); j++)

    // Scatter contact forces to the nodes
    // WARNING: the following scatter operation will not yield parallel consistency with multi-threading!!!
    for (int i=0; i<Nnodes; i++) {
      f[Ndofs_per_node*nodes[i]+0] -= f_nodes[Ndofs_per_node*i+0];
      f[Ndofs_per_node*nodes[i]+1] -= f_nodes[Ndofs_per_node*i+1];
    }
    for (int i=0; i<Nsegment_nodes; i++) {
      f[Ndofs_per_node*segment_nodes[i]+0] -= f_segment_nodes[Ndofs_per_node*i+0];
      f[Ndofs_per_node*segment_nodes[i]+1] -= f_segment_nodes[Ndofs_per_node*i+1];
    }
    
  } // update_contact_forces

}; // ContactInteraction

#endif // CONTACT_INTERACTION_H
