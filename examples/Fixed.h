#ifndef FIXED_H
#define FIXED_H

#include "types.h"
#include <iostream>
#include <math.h>

template<int RADIX = 10, int EXPONENT = 0>
struct Fixed {
  inline static double decimal = std::pow(RADIX,EXPONENT);
  Integer mantissa;

  // constructor for a Fixed-precision number given its mantissa
  Fixed(Integer i = 0) : mantissa(i) { }

  // constructor for a Fixed-precision number given a floating point number
  Fixed(double d = 0.0) : mantissa(d/decimal) { }

  // implicit type-cast to a floating point number
  operator double() { return mantissa*decimal; }

  // summation of two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator+(Fixed<RADIX,EXPONENT> const summand) {
    return Fixed<RADIX,EXPONENT>(mantissa + summand.mantissa);
  } // operator+

  // subtraction of two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator-(Fixed<RADIX,EXPONENT> const subtrahend) {
    return Fixed<RADIX,EXPONENT>(mantissa - subtrahend.mantissa);
  } // operator-

  // multiplication of a Fixed-precision number by an Integer
  Fixed<RADIX,EXPONENT> operator*(Integer const multiplier) {
    return Fixed<RADIX,EXPONENT>(mantissa*multiplier);
  } // operator*

  // division of a Fixed-precision number by an Integer
  Fixed<RADIX,EXPONENT> operator/(Integer const divisor) {
    return Fixed<RADIX,EXPONENT>(mantissa/divisor);
  } // operator/

  // modulo (remainder) of a Fixed-precision number
  Fixed<RADIX,EXPONENT> operator%(Integer const divisor) {
    return Fixed<RADIX,EXPONENT>(mantissa%divisor);
  } // operator%
  
}; // Fixed

#endif // FIXED_H
