#include "gtest/gtest.h"
#include "Rational.h"
#include <iostream>

TEST(test_Rational, constructors) {
  // Test creation of Rational numbers
  Rational a = Rational(  5,6);
  Rational b = Rational(-17,12);
  ASSERT_EQ(a, +5.0/6.0);
  ASSERT_EQ(b,-17.0/12.0);
} /* TEST(test_Rational, constructors) */

TEST(test_Rational, arithmetic) {
  Rational a = Rational(  5,6);
  Rational b = Rational(-17,12);
  Integer  i = 5;
  // Test basic arithmetic operations between two Rational numbers
  ASSERT_EQ(a+b, -7.0/12.0);
  ASSERT_EQ(a-b,+27.0/12.0);
  // Test basic arithemtic operations between a Rational number and an Integer
  ASSERT_EQ(a*i,  (5*i)/6.0 );
  ASSERT_EQ(b*i,(-17*i)/12.0);
  ASSERT_EQ(a/i,  (5/i)/6.0 );
  ASSERT_EQ(b/i,(-17/i)/12.0);
  ASSERT_EQ(a%i,  (5%i)/6.0 );
  ASSERT_EQ(b%i,(-17%i)/12.0);
} /* TEST(test_Rational, arithmetic) */

TEST(test_Rational, comparison) {
  // Test comparison of Rational numbers
  Rational a = Rational(  5,6);
  Rational b = Rational(-17,12);
  ASSERT_TRUE (a >  b);
  ASSERT_FALSE(a <  b);
  ASSERT_TRUE (a >= b);
  ASSERT_FALSE(a <= b);
  ASSERT_FALSE(a == b);
  ASSERT_TRUE (a != b);
} /* TEST(test_Rational, comparison) */

TEST(test_Rational, cout) {
  // Test printing of Rational numbers
  Rational a = Rational(  5,6);
  Rational b = Rational(-17,12);
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "a + b = " << (a+b) << std::endl;
  std::cout << "a - b = " << (a-b) << std::endl;
} /* TEST(test_Rational, cout) */
