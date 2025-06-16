#ifndef MATRIX_H
#define MATRIX_H

#include "types.h"
#include <vector>
#include <iostream>
#include <math.h>
#include "Dual.h"

// Symmetrix NxN matrix (assumed positive-definite)
template<int N>
struct Matrix {
  double data[N][N];

  Matrix(double initializer = 0.0) {
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        data[i][j] = initializer;
      }
    }
  }

  Matrix(Matrix<N>& const copy) {
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        data[i][j] = copy(i,j);
      }
    }
  }

  // Zero all entries in the matrix
  void zero(void) {
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	data[i][j] = 0.0;
      }
    }
  } // zero()

  // Perform LDL factorization
  void factor(Matrix<N>& L, double* D) const {
    // Loop over all j rows
    for (int j=0; j<N; j++) {
      D[j] = data[j][j];
      // Loop over all k columns preceding column j
      for (int k=0; k<j; k++) {
	D[j] -= L(j,k)*L(j,k)*D[k];
      }
      L(j,j) = 1.0;
      // Loop over all i rows below row j
      for (int i=j+1; i<N; i++) {
	L(i,j) = data[i][j];
	// Loop over all k columns preceding column j
	for (int k=0; k<j; k++) {
	  L(i,j) -= L(i,k)*L(j,k)*D[k];
	}
	L(i,j) /= D[j];
	L(j,i) = 0.0;
      }
    }
  } // factor()

  // Transpose matrix
  void transpose(void) {
    for (int i=0; i<N; i++) {
      for (int j=0; j<i; j++) {
	double aij = data[i][j];
	data[i][j] = data[j][i];
	data[j][i] = aij;
      }
    }
  } // transpose()

  template<class T>
  void lower_dual_update(std::vector<Dual<T> >& const x, std::vector<Dual<T> >& y) const {
    for (int i=0; i<N; i++) {
      y[i].first = x[i].first;
      for (int j=0; j<i; j++) {
        y[i].first += data[i][j]*x[j].first;
      }
    }
    for (int i=N-1; i>=0; i--) {
      y[i].second = x[i].second;
      for (int j=i+1; j<N; j++) {
        y[i].second -= data[j][i]*y[j].second;
      }
    }
  } // lower_dual_update

  template<class T>
  void upper_dual_update(std::vector<Dual<T> >& const x, std::vector<Dual<T> >& y) const {
    for (int i=N-1; i>=0; i--) {
      y[i].first = x[i].first;
      for (int j=i+1; j<N; j++) {
        y[i].first += data[i][j]*x[j].first;
      }
    }
    for (int i=0; i<N; i++) {
      y[i].second = x[i].second;
      for (int j=0; j<i; j++) {
        y[i].second -= data[j][i]*y[j].second;
      }
    }
  } // upper_dual_update

  // Access operator
  double& operator()(int i, int j)       { return data[i][j]; }
  double  operator()(int i, int j) const { return data[i][j]; }

}; // Matrix

template<class T, int N>
void dual_update(Matrix<N>& L, Matrix<N>& U, double* D, std::vector<Dual<T> >& const x, std::vector<Dual<T> >& y) {
  std::vector<Dual<T> > temp(N);
  U.upper_dual_update(x,temp);
  for (int i=0; i<N; i++) temp[i] = D[i] * temp[i];
  L.lower_dual_update(temp,y);
} // dual_update

// overloaded output stream operator
template<int N>
std::ostream& operator<<(std::ostream& os, const Matrix<N>& mat) {
  os << "[ ";
  for (int i=0; i<N; i++) {
    if (i > 0) os << std::endl;
    os << "[ ";
    for (int j=0; j<N; j++) {
      os << mat(i,j) << ' ';
    }
    os << "]";
  }
  os << " ]";
  return os;
} // operator<<

#endif // MATRIX_H
