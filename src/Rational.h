#ifndef RATIONAL_H
#define RATIONAL_H

#include "types.h"
#include <iostream>

struct Rational {
  Integer numerator;
  Integer denominator;

  // basic constructor for a Rational number given two Integer integer numbers
  Rational(Integer n = 0, Integer d = 1) : numerator(n), denominator(d) { }

  // basic constructor for a Rational number given two Integer integer numbers
  //Rational(double f, Integer d) : numerator(f*d), denominator(d) { }

  // implicit type-cast to a floating point number
  operator double() { return double(numerator)/denominator; }

  // summation of two Rational numbers
  Rational operator+(Rational const summand) {
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
  Rational operator-(Rational const subtrahend) {
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
  Rational operator*(Integer const multiplier) { return Rational(numerator*multiplier,denominator); }
  Rational operator/(Integer const divisor)    { return Rational(numerator/divisor,   denominator); }
  Rational operator%(Integer const divisor)    { return Rational(numerator%divisor,   denominator); }

  // comparison of two Fixed-precision numbers
  bool operator> (Rational const other) { return (this->double() >  double(other)); }
  bool operator< (Rational const other) { return (this->double() <  double(other)); }
  bool operator>=(Rational const other) { return (this->double() >= double(other)); }
  bool operator<=(Rational const other) { return (this->double() <= double(other)); }
  bool operator==(Rational const other) { return (this->double() == double(other)); }
  bool operator!=(Rational const other) { return (this->double() != double(other)); }
  
}; // Rational

// overloaded output stream operator
std::ostream& operator<<(std::ostream& os, const Rational& number) {
  os << number.numerator << '/' << number.denominator;
  return os;
} // operator<<

#endif // RATIONAL_H
