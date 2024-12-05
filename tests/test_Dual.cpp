#include "gtest/gtest.h"
#include "Fixed.h"
#include "Rational.h"
#include "Dual.h"
#include <iostream>

// Declare a standard Fixed-precision number and correponding Dual number type
const int          RADIX = 10;
const int       EXPONENT = -6;
const double     decimal = std::pow(RADIX,+EXPONENT);
const int    inv_decimal = std::pow(RADIX,-EXPONENT);
typedef Fixed<RADIX,EXPONENT> FixedE6;
typedef Dual<FixedE6>         DualE6;

TEST(test_Dual, constructors) {
  // Test creation of Dual numbers
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  DualE6  d = DualE6(a,b);
  ASSERT_EQ(d.first, a);
  ASSERT_EQ(d.second,b);
} /* TEST(test_Dual, constructors) */

TEST(test_Dual, arithmetic) {
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  DualE6  c = DualE6(-b,a);
  DualE6  d = DualE6(a,b);
  Integer i = 5;
  // Test basic arithmetic operations between two Dual numbers
  ASSERT_EQ(c+d,DualE6(-b+a,a+b));
  ASSERT_EQ(c-d,DualE6(-b-a,a-b));
  // Test basic arithemtic operations between a Dual number and an Integer
  ASSERT_EQ(c*i,DualE6(-b*i+a%i,a/i));
  ASSERT_EQ(c/i,DualE6(-b/i,a*i-b%i));
  ASSERT_EQ(d*i,DualE6(a*i-b%i,b/i));
  ASSERT_EQ(d/i,DualE6(a/i,b*i-a%i));
} /* TEST(test_Dual, arithmetic) */

TEST(test_Dual, reversibility) {
  FixedE6  a = +5.0/6.0;
  FixedE6  b =-17.0/12.0;
  DualE6   c = DualE6(-b,a);
  DualE6   d = DualE6(a,b);
  Integer  i = 5;
  Rational r = Rational(+5,6);
  // Test reversibility of basic arithmetic operations between two Dual numbers
  ASSERT_EQ((c+d)-c,d);
  ASSERT_EQ((c-d)+d,c);
  // Test reversibility of basic arithemtic operations between a Dual number and an Integer
  ASSERT_EQ((c/i)*i,c);
  ASSERT_EQ((c*i)/i,c);
  ASSERT_EQ((d/i)*i,d);
  ASSERT_EQ((d*i)/i,d);
  // Test reversibility of basic arithmetic operations between a Dual number and a Rational number
  ASSERT_EQ((c/r)*r,c);
  ASSERT_EQ((c*r)/r,c);
  ASSERT_EQ((d/r)*r,d);
  ASSERT_EQ((d*r)/r,d);
} /* TEST(test_Dual, reversibility) */

TEST(test_Dual, comparison) {
  // Test comparison of Dual numbers
  FixedE6 a = +5.0/6.0;
  FixedE6 b =-17.0/12.0;
  DualE6  c = DualE6(b,a);
  DualE6  d = DualE6(a,b);
  ASSERT_TRUE (c == c);
  ASSERT_FALSE(c == d);
  ASSERT_FALSE(c != c);
  ASSERT_TRUE (c != d);
} /* TEST(test_Dual, comparison) */

TEST(test_Dual, cout) {
  // Test printing of Dual numbers
  FixedE6  a = +5.0/6.0;
  FixedE6  b =-17.0/12.0;
  DualE6   c = DualE6(-b,a);
  DualE6   d = DualE6(a,b);
  Integer  i = 5;
  Rational r = Rational(+5,6);
  std::cout << "c = " << c << std::endl;
  std::cout << "d = " << d << std::endl;
  std::cout << "c + d = " << (c+d) << std::endl;
  std::cout << "c - d = " << (c-d) << std::endl;
  std::cout << "(c / i) * i = " << (c/i)*i << std::endl;
  std::cout << "(c * i) / i = " << (c*i)/i << std::endl;
  std::cout << "(d / i) * i = " << (d/i)*i << std::endl;
  std::cout << "(d * i) / i = " << (d*i)/i << std::endl;
  std::cout << "(c / r) * r = " << (c/r)*r << std::endl;
  std::cout << "(c * r) / r = " << (c*r)/r << std::endl;
  std::cout << "(d / r) * r = " << (d/r)*r << std::endl;
  std::cout << "(d * r) / r = " << (d*r)/r << std::endl;
} /* TEST(test_Dual, cout) */
