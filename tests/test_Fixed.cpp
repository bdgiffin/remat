#include "gtest/gtest.h"
#include "Fixed.h"
#include <iostream>
#include <math.h>

// Declare a standard Fixed-precision number
const int          RADIX = 10;
const int       EXPONENT = -6;
const double     decimal = std::pow(RADIX,+EXPONENT);
const int    inv_decimal = std::pow(RADIX,-EXPONENT);
typedef Fixed<RADIX,EXPONENT> FixedE6;

TEST(test_Fixed, constructors) {
  // Test creation of Fixed-precision numbers
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  ASSERT_EQ(a.mantissa, +833333);
  ASSERT_EQ(b.mantissa,-1416666);
} /* TEST(test_Fixed, constructors) */

TEST(test_Fixed, arithmetic) {
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  Integer i = 5;
  // Test basic arithmetic operations between two Fixed-precision numbers
  ASSERT_EQ((a+b).mantissa,- 583333);
  ASSERT_EQ((a-b).mantissa,+2249999);
  // Test basic arithemtic operations between a Fixed-precision number and an Integer
  ASSERT_EQ((a*i).mantissa,+4166665);
  ASSERT_EQ((b*i).mantissa,-7083330);
  ASSERT_EQ((a/i).mantissa,+ 166666);
  ASSERT_EQ((b/i).mantissa,- 283334);
  ASSERT_EQ((a%i).mantissa,+      3);
  ASSERT_EQ((b%i).mantissa,+      4);
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

