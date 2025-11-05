#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <vector>
#include "types.h"

// Collection of node identifiers that share a common boundary condition.
struct NodeSet {
  std::vector<int> node_ids;
};

// Define function pointer for time-varying displacement boundary conditions.
typedef Real (*TimeFunction)(Real time, Real x, Real y, int node_id, int component);

// Data associated with a prescribed, time-varying displacement boundary condition.
struct DisplacementBoundaryCondition {
  NodeSet nodes;
  int component = 0;
  TimeFunction function = nullptr;
  std::vector<Real> last_values;
  std::vector<Real> prescribed_velocities;
};

#endif // BOUNDARY_CONDITION_H
