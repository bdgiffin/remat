#ifndef FIXED_H
#define FIXED_H

#include "types.h"
#include <iostream>
#include <math.h>

template<int RADIX = 10, int EXPONENT = 0>
struct Fixed {
  inline static double decimal = std::pow(RADIX,EXPONENT);
  Integer mantissa;

  // constructor methods for a Fixed-precision number
  Fixed(double d = 0.0) : mantissa(d/decimal) { }

  // implicit type-cast to a floating point number
  operator double() { return mantissa*decimal; }

  // basic arithmetic operations between two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator+(Fixed<RADIX,EXPONENT> const summand)    { return Fixed<RADIX,EXPONENT>(mantissa + summand.mantissa);    }
  Fixed<RADIX,EXPONENT> operator-(Fixed<RADIX,EXPONENT> const subtrahend) { return Fixed<RADIX,EXPONENT>(mantissa - subtrahend.mantissa); }

  // basic arithmetic operations between a Fixed-precision number and an Integer
  Fixed<RADIX,EXPONENT> operator*(Integer const multiplier) { return Fixed<RADIX,EXPONENT>(mantissa*multiplier); }
  Fixed<RADIX,EXPONENT> operator/(Integer const divisor)    { return Fixed<RADIX,EXPONENT>(mantissa/divisor);    }
  Fixed<RADIX,EXPONENT> operator%(Integer const divisor)    { return Fixed<RADIX,EXPONENT>(mantissa%divisor);    }

  // comparison of two Fixed-precision numbers
  bool operator> (Fixed<RADIX,EXPONENT> const other) { return (mantissa >  other.mantissa); }
  bool operator< (Fixed<RADIX,EXPONENT> const other) { return (mantissa <  other.mantissa); }
  bool operator>=(Fixed<RADIX,EXPONENT> const other) { return (mantissa >= other.mantissa); }
  bool operator<=(Fixed<RADIX,EXPONENT> const other) { return (mantissa <= other.mantissa); }
  bool operator==(Fixed<RADIX,EXPONENT> const other) { return (mantissa == other.mantissa); }
  bool operator!=(Fixed<RADIX,EXPONENT> const other) { return (mantissa != other.mantissa); }
  
}; // Fixed

#endif // FIXED_H
