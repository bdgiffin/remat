#ifndef ELEMENT_H
#define ELEMENT_H

#include "Parameters.h"
#include <vector>
#include <math.h>

template<typename Material_T>
class Element {
public:
  Material_T m_model; // material model

  // Empty constructor
  Element(void) { }

  // Parameterized constructor
  Element(Parameters& params) : m_model(params) { }
    
  // Return the number of state variables for allocation purposes
  int num_state_vars(void) { return (4+1)*m_model.num_state_vars(); }

  // Initialize the element state
  void initialize(Real* state) {
    
    // zero-initialize element-averaged state
    // (stored at the fictitious integration point 0)
    const int num_state_vars = m_model.num_state_vars();
    for (int i=0; i<num_state_vars; i++) {
      state[i] = 0.0;
    }
    
    // loop over integration points
    // (stored at the actual integration point 1-4)
    for (int q=0; q < 4; q++) {
      Real* model_state = &state[(q+1)*num_state_vars];
      m_model.initialize(model_state);

      // sum contributions to element-averaged state
      for (int i=0; i<num_state_vars; i++) {
	state[i] += 0.25*model_state[i];
      }
    }
    
  } // initialize()
    
  // Initialize variable material properties for the current element
  void initialize_variable_properties(Real (&x)[8], Real* state, double (*function_xy)(double,double)) {
    
    const Real sqrt_third = 1.0/std::sqrt(3.0);
      
    // Define parent nodal coordinates
    Real xi[4]   = { -1.0, +1.0, +1.0, -1.0 };
    Real eta[4]  = { -1.0, -1.0, +1.0, +1.0 };

    // zero-initialize element-averaged state
    const int num_state_vars = m_model.num_state_vars();
    for (int i=0; i<num_state_vars; i++) {
      state[i] = 0.0;
    }
      
    // loop over integration points
    for (int q = 0; q < 4; q++) {

      // Set quadrature point position
      Real xiq  = xi[q]*sqrt_third;
      Real etaq = eta[q]*sqrt_third;
	
      // Compute shape function values at the current quadrature point (phi)
      Real phi[4];
      for (int i = 0; i < 4; i++) {
	phi[i] = 0.25*(1.0+xiq*xi[i])*(1.0+etaq*eta[i]);
      }

      // Compute interpolated coordinates of the current quadrature point (xq)
      Real xq[2] = {0.0, 0.0};
      for (int i = 0; i < 4; i++) {
	xq[0] += x[2*i+0] * phi[i];
	xq[1] += x[2*i+1] * phi[i];
      }

      // Initialize variable material properties for the current material point
      Real* model_state = &state[(q+1)*num_state_vars];
      Real psi;
      m_model.initialize_variable_properties(xq,model_state,function_xy);

      // sum contributions to element-averaged state
      for (int i=0; i<num_state_vars; i++) {
	state[i] += 0.25*model_state[i];
      }
	
    }
  } // initialize_variable_properties()
    
  // Update the element state using the current nodal displacements
  void update(Real (&x)[8], Real (&u)[8], Real (&m)[8], Real (&f)[8], Real &E, Real* state, Real dt) {

    // Determine nodal masses:
    {
      // Compute mid-point shape function partial derivatives (dphi/dxi = 0.25*xi)
      Real dxi[4]  = { -0.25, +0.25, +0.25, -0.25 };
      Real deta[4] = { -0.25, -0.25, +0.25, +0.25 };

      // Compute mid-point Jacobian (J = dx/dxi = sum_a x_a (x) dphi/dxi)
      Real J[2][2] = { {0.0,0.0}, {0.0,0.0} };
      for (int i = 0; i < num_nodes(); i++) {
	J[0][0] += x[2*i+0] * dxi[i]; J[0][1] += x[2*i+0] * deta[i];
	J[1][0] += x[2*i+1] * dxi[i]; J[1][1] += x[2*i+1] * deta[i];
      }

      // Compute determinant of J
      Real detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

      // Compute total element mass going to each node = (4 x detJ x density) / 4
      Real mass = detJ*m_model.density();

      // Sum nodal mass contributions
      for (int i = 0; i < 2*num_nodes(); i++) m[i] = mass;
    }

    // Determine nodal forces:
    {
      const Real sqrt_third = 1.0/std::sqrt(3.0);
      
      // Define parent nodal coordinates
      Real xi[4]   = { -1.0, +1.0, +1.0, -1.0 };
      Real eta[4]  = { -1.0, -1.0, +1.0, +1.0 };

      // zero-initialize element-averaged state
      const int num_state_vars = m_model.num_state_vars();
      for (int i=0; i<num_state_vars; i++) {
	state[i] = 0.0;
      }
      
      // loop over integration points
      for (int q = 0; q < 4; q++) {

	// Set quadrature point position
	Real xiq   = xi[q]*sqrt_third;
	Real etaq  = eta[q]*sqrt_third;
	
	// Compute shape function partial derivatives (dphi/dxi)
	Real dxi[4], deta[4];
	for (int i = 0; i < 4; i++) {
	  dxi[i]  = 0.25*xi[i]*(1.0+etaq*eta[i]);
	  deta[i] = 0.25*eta[i]*(1.0+xiq*xi[i]);
	}

	// Compute initial and current Jacobians (J = dx/dxi = sum_a x_a (x) dphi/dxi)
	Real J0[2][2] = { {0.0,0.0}, {0.0,0.0} };
	Real J[2][2]  = { {0.0,0.0}, {0.0,0.0} };
	for (int i = 0; i < 4; i++) {
	  J0[0][0] += x[2*i+0] * dxi[i]; J0[0][1] += x[2*i+0] * deta[i];
	  J0[1][0] += x[2*i+1] * dxi[i]; J0[1][1] += x[2*i+1] * deta[i];
	  Real xt = x[2*i+0] + u[2*i+0]; Real yt = x[2*i+1] + u[2*i+1];
	  J[0][0] += xt * dxi[i]; J[0][1] += xt * deta[i];
	  J[1][0] += yt * dxi[i]; J[1][1] += yt * deta[i];
	}

	// Compute determinant of J0 and its inverse
	Real detJ0 = J0[0][0]*J0[1][1] - J0[0][1]*J0[1][0];
	Real inv_detJ0 = 1.0/detJ0;

	// Compute inv(J0)
	Real invJ0[2][2];
	invJ0[0][0] = +J0[1][1]*inv_detJ0;
	invJ0[0][1] = -J0[0][1]*inv_detJ0;
	invJ0[1][0] = -J0[1][0]*inv_detJ0;
	invJ0[1][1] = +J0[0][0]*inv_detJ0;
	
	// Compute the current deformation gradient F = J*invJ0
	Real F[2][2];
	F[0][0] = J[0][0]*invJ0[0][0] + J[0][1]*invJ0[1][0];
	F[0][1] = J[0][0]*invJ0[0][1] + J[0][1]*invJ0[1][1];
	F[1][0] = J[1][0]*invJ0[0][0] + J[1][1]*invJ0[1][0];
	F[1][1] = J[1][0]*invJ0[0][1] + J[1][1]*invJ0[1][1];

	// update the material state
	Real* model_state = &state[(q+1)*num_state_vars];
	Real psi;
	m_model.update(F,psi,model_state,dt);

	// sum contributions to element-averaged state
	for (int i=0; i<num_state_vars; i++) {
	  state[i] += 0.25*model_state[i];
	}

	// Compute cof(J)
	Real cofJ[2][2];
	cofJ[0][0] = +J[1][1];
	cofJ[1][0] = -J[0][1];
	cofJ[0][1] = -J[1][0];
	cofJ[1][1] = +J[0][0];

	// Compute P = sigma*cofJ relative to the parent element configuration
	Real P[2][2];
	P[0][0] = model_state[0]*cofJ[0][0] + model_state[5]*cofJ[1][0];
	P[0][1] = model_state[0]*cofJ[0][1] + model_state[5]*cofJ[1][1];
	P[1][0] = model_state[5]*cofJ[0][0] + model_state[1]*cofJ[1][0];
	P[1][1] = model_state[5]*cofJ[0][1] + model_state[1]*cofJ[1][1];

	// sum nodal forces by stress divergence: P.dphi/dxi d[]
	for (int i = 0; i < num_nodes(); i++) {
	  f[2*i+0] += P[0][0]*dxi[i] + P[0][1]*deta[i];
	  f[2*i+1] += P[1][0]*dxi[i] + P[1][1]*deta[i];
	}

	// sum contribution to the total elastic strain energy
	E += psi*detJ0;
	
      } // end loop over integration points
    }

  } // update()

  // Conditionally load material history parameters from memory for this element
  void load_state(Real* state, std::vector<Real>& overflow_state) {
    const int num_state_vars = m_model.num_state_vars();

    // loop over integration points (in reverse order)
    for (int q=3; q>=0; q--) {
      m_model.load_state(&state[num_state_vars*(q+1)],overflow_state);
    }
  }

  // Conditionally store material history parameters in memory for this element
  void store_state(Real* state, std::vector<Real>& overflow_state) {
    const int num_state_vars = m_model.num_state_vars();

    // loop over integration points (in forward order)
    for (int q=0; q<4; q++) {
      m_model.store_state(&state[num_state_vars*(q+1)],overflow_state);
    }
  }
    
  // Return the number of nodes per element
  constexpr int num_nodes(void) { return 4; }

}; // Element

#endif // ELEMENT_H
