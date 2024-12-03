#include "gtest/gtest.h"
#include "Fixed.h"
#include <iostream>
#include <math.h>

// Declare a standard Fixed-precision number
const int RADIX = 10;
const int EXPONENT = -6;
const double decimal = std::pow(RADIX,EXPONENT)
typedef Fixed<RADIX,EXPONENT> FixedE6;

TEST(test_Fixed, constructors) {
  // Test creation of Fixed-precision numbers
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  ASSERT_EQ(a,+0.833333);
  ASSERT_EQ(b,-1.416667);
} /* TEST(test_Fixed, constructors) */

TEST(test_Fixed, arithmetic) {
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  Integer i = 5;
  // Test basic arithmetic operations between two Fixed-precision numbers
  ASSERT_EQ(a+b,-0.583334);
  ASSERT_EQ(a-b,+2.250000);
  // Test basic arithemtic operations between a Fixed-precision number and an Integer
  ASSERT_EQ(a*i,+4.166665);
  ASSERT_EQ(b*i,-7.083335);
  ASSERT_EQ(a/i,+0.166666);
  ASSERT_EQ(b/i,-0.283333);
  ASSERT_EQ(a%i,+0.000003);
  ASSERT_EQ(b%i,+0.000003);
} /* TEST(test_Fixed, arithmetic) */

TEST(test_Fixed, comparison) {
  // Test creation of Fixed-precision numbers
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  // Test creation of Fixed-precision numbers
  ASSERT_TRUE (a >  b);
  ASSERT_FALSE(a <  b);
  ASSERT_TRUE (a >= b);
  ASSERT_FALSE(a <= b);
  ASSERT_FALSE(a == b);
  ASSERT_TRUE (a != b);
} /* TEST(test_Fixed, comparison) */

