#ifndef MATRIX_H
#define MATRIX_H

#include "types.h"
#include "Rational.h"
#include <iostream>
#include <math.h>

// Symmetrix NxN matrix (assumed positive-definite)
template<int N>
struct Matrix {
  double A[N][N];

  Matrix(void) {}

  void factor(void) {
    double L[N][N];
    double D[N];
    // Loop over all j rows
    for (int j=0; j<N; j++) {
      D[j] = A[j][j];
      // Loop over all k columns preceding column j
      for (int k=0; k<j; k++) {
	D[j] -= L[j][k]*L[j][k]*D[k];
      }
      // Loop over all i rows below row j
      for (int i=j+1; i<N; i++) {
	L[i][j] = A[i][j];
	// Loop over all k columns preceding column j
	for (int k=0; k<j; k++) {
	  L[i][j] -= L[i][k]*L[j][k]*D[k];
	}
	L[i][j] /= D[j];
      }
    }
  } // factor()

}; // Matrix
