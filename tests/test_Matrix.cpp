#include "gtest/gtest.h"
#include "Matrix.h"
#include <iostream>

TEST(test_Matrix, constructors) {
  // Test creation of matrices
  const int N = 10;
  Matrix<N> A;
  Matrix<N> B(2.0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      ASSERT_EQ(A(i,j), 0.0);
      ASSERT_EQ(B(i,j), 2.0);
    }
  }
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "B = " << std::endl << B << std::endl;
} /* TEST(test_Matrix, constructors) */

TEST(test_Matrix, factor) {
  // Test simple LDL factorization
  const int N = 10;
  Matrix<N> A,L;
  double D[N];
  for (int i=0; i<N; i++) {
    A(i,i) = 2.0;
    if (i > 0)   A(i,i-1) = -1.0;
    if (i < N-1) A(i,i+1) = -1.0;
  }
  A.factor(L,D);

  // Compare against reference solution obtained from MATLAB
  double tolerance = 1.0e-4;
  double Lref[N] = {-0.5000,-0.6667,-0.7500,-0.8000,-0.8333,-0.8571,-0.8750,-0.8889,-0.9000,0.0000};
  double Dref[N] = {2.0000,1.5000,1.3333,1.2500,1.2000,1.1667,1.1429,1.1250,1.1111,1.1000};
  for (int i=0; i<N; i++) {
    ASSERT_NEAR(D[i], Dref[i], tolerance);
    ASSERT_EQ(L(i,i), 1.0);
    for (int j=0; j<i; j++) ASSERT_EQ(L(j,i), 0.0);
    if (i < (N-1)) ASSERT_NEAR(L(i+1,i), Lref[i], tolerance);
    for (int j=i+2; j<N; j++) ASSERT_EQ(L(j,i), 0.0);
  }

  // Optionally print out the results
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "L = " << std::endl << L << std::endl;
  std::cout << "D = [ ";
  for (int i=0; i<N; i++) {
    std::cout << D[i] << " ";
  }
  std::cout << "]" << std::endl;
  
} /* TEST(test_Matrix, factor) */

TEST(test_Matrix, transpose) {
  // LDL factorization
  const int N = 10;
  Matrix<N> A,L,U;
  double D[N];
  for (int i=0; i<N; i++) {
    A(i,i) = 2.0;
    if (i > 0)   A(i,i-1) = -1.0;
    if (i < N-1) A(i,i+1) = -1.0;
  }
  A.factor(L,D);

  // Test transpose
  U = L;
  U.transpose;

  // Compare against reference solution obtained from MATLAB
  double tolerance = 1.0e-4;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) ASSERT_EQ(L(j,i), U(i,j));
  }

  // Optionally print out the results
  std::cout << "L = " << std::endl << L << std::endl;
  std::cout << "U = " << std::endl << U << std::endl;
  
} /* TEST(test_Matrix, transpose) */

TEST(test_Matrix, transpose) {
  // LDL factorization
  const int N = 10;
  Matrix<N> A,L,U;
  double D[N];
  for (int i=0; i<N; i++) {
    A(i,i) = 2.0;
    if (i > 0)   A(i,i-1) = -1.0;
    if (i < N-1) A(i,i+1) = -1.0;
  }
  A.factor(L,D);

  // Test transpose
  U = L;
  U.transpose;

  // Compare against reference solution obtained from MATLAB
  double tolerance = 1.0e-4;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) ASSERT_EQ(L(j,i), U(i,j));
  }

  // Optionally print out the results
  std::cout << "L = " << std::endl << L << std::endl;
  std::cout << "U = " << std::endl << U << std::endl;
  
} /* TEST(test_Matrix, transpose) */
