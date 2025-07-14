#ifndef RATIONAL_H
#define RATIONAL_H

#include "types.h"
#include <iostream>
#include <math.h>
#include <limits>

struct Rational {
  Integer numerator;
  Integer denominator;

  // default constructor
  Rational(void) : numerator(0), denominator(1) { }

  // basic constructor for a Rational number given two Integer integer numbers
  Rational(Integer n, Integer d) : numerator(n), denominator(d) { }

  // basic constructor for a Rational number given a Real (floating point) number
  Rational(Real f) {
    //// Separate exponent and mantissa in [.5, 1),
    //// and scale the mantissa by the number of digits to produce the integer numerator
    //int exponent;
    //numerator = scalb(std::frexp(f, &exponent), std::numeric_limits<Real>::digits);
    //
    //// Adjust the exponent to compensate for scaling of the mantissa
    //exponent -= std::numeric_limits<Real>::digits;
    //
    //// Scale either the numerator or denominator, depending on the sign of the exponent
    //denominator = 1;
    //if (exponent < 0) {
    //  denominator = scalb(denominator, -exponent);
    //} else {
    //  numerator = scalb(numerator, exponent);
    //}

    // Choose a standard denominator relative to the largest Integer
    denominator = std::sqrt(std::numeric_limits<Integer>::max());
    numerator   = Integer(f*denominator);
  }

  // implicit type-cast to a floating point number
  Real decimal(void) const { return Real(numerator)/denominator; }
  operator Real()    const { return decimal(); }

  // summation of two Rational numbers
  Rational operator+(Rational const summand) const {
    Rational sum;
    if (denominator == summand.denominator) {
      sum.numerator   = numerator + summand.numerator;
      sum.denominator = denominator;
    } else if ((denominator > summand.denominator) && (denominator % summand.denominator == 0)) {
      sum.numerator   = numerator + summand.numerator * (denominator / summand.denominator);
      sum.denominator = denominator;
    } else if ((summand.denominator > denominator) && (summand.denominator % denominator == 0)) {
      sum.numerator   = numerator * (summand.denominator / denominator) + summand.numerator;
      sum.denominator = summand.denominator;
    } else {
      sum.numerator   = numerator * summand.denominator + summand.numerator * denominator;
      sum.denominator = denominator * summand.denominator;
    }
    return sum;
  } // operator+

  // subtraction of two Rational numbers
  Rational operator-(Rational const subtrahend) const {
    Rational difference;
    if (denominator == subtrahend.denominator) {
      difference.numerator   = numerator - subtrahend.numerator;
      difference.denominator = denominator;
    } else if ((denominator > subtrahend.denominator) && (denominator % subtrahend.denominator == 0)) {
      difference.numerator   = numerator - subtrahend.numerator * (denominator / subtrahend.denominator);
      difference.denominator = denominator;
    } else if ((subtrahend.denominator > denominator) && (subtrahend.denominator % denominator == 0)) {
      difference.numerator   = numerator * (subtrahend.denominator / denominator) - subtrahend.numerator;
      difference.denominator = subtrahend.denominator;
    } else {
      difference.numerator   = numerator * subtrahend.denominator - subtrahend.numerator * denominator;
      difference.denominator = denominator * subtrahend.denominator;
    }
    return difference;
  } // operator-
  
  // basic arithmetic operations between a Rational number and an Integer
  Rational operator*(Integer const multiplier) const { return Rational(numerator*multiplier,denominator); }
  Rational operator/(Integer const divisor)    const { return Rational(numerator/divisor,   denominator); }
  Rational operator%(Integer const divisor)    const { return Rational(numerator%divisor,   denominator); }

  // comparison of two Fixed-precision numbers
  bool operator> (Rational const other) const { return (decimal() >  other.decimal()); }
  bool operator< (Rational const other) const { return (decimal() <  other.decimal()); }
  bool operator>=(Rational const other) const { return (decimal() >= other.decimal()); }
  bool operator<=(Rational const other) const { return (decimal() <= other.decimal()); }
  bool operator==(Rational const other) const { return (decimal() == other.decimal()); }
  bool operator!=(Rational const other) const { return (decimal() != other.decimal()); }
  
}; // Rational

// overloaded output stream operator
std::ostream& operator<<(std::ostream& os, const Rational& number) {
  os << number.numerator << '/' << number.denominator;
  return os;
} // operator<<

#endif // RATIONAL_H
