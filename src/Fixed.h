#ifndef FIXED_H
#define FIXED_H

#include "types.h"
#include <iostream>
#include <math.h>

template<int RADIX = 10, int EXPONENT = 0>
struct Fixed {
  inline static double     decimal = std::pow(RADIX,+EXPONENT);
  inline static double inv_decimal = std::pow(RADIX,-EXPONENT);
  Integer mantissa;

  // constructor method for a Fixed-precision number
  Fixed(double d = 0.0) : mantissa(d*inv_decimal) { }

  // implicit type-cast to a floating point number
  operator double() const { return mantissa*decimal; }

  // basic arithmetic operations between two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator+(Fixed<RADIX,EXPONENT> const summand)    const { return Fixed<RADIX,EXPONENT>(mantissa + summand.mantissa);    }
  Fixed<RADIX,EXPONENT> operator-(Fixed<RADIX,EXPONENT> const subtrahend) const { return Fixed<RADIX,EXPONENT>(mantissa - subtrahend.mantissa); }

  // basic arithmetic operations between a Fixed-precision number and an Integer
  Fixed<RADIX,EXPONENT> operator*(Integer const multiplier) const { return Fixed<RADIX,EXPONENT>(mantissa*multiplier); }
  Fixed<RADIX,EXPONENT> operator/(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(mantissa/divisor);    }
  Fixed<RADIX,EXPONENT> operator%(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(mantissa%divisor);    }

  // basic arithmetic operations between a Fixed-precision number and a double-precision floating point number
  Fixed<RADIX,EXPONENT> operator*(double const multiplier) const { return Fixed<RADIX,EXPONENT>(Integer(mantissa*multiplier)); }
  Fixed<RADIX,EXPONENT> operator/(double const divisor)    const { return Fixed<RADIX,EXPONENT>(Integer(mantissa/divisor));    }
  Fixed<RADIX,EXPONENT> operator%(double const divisor)    const { return Fixed<RADIX,EXPONENT>(Integer(std::remainder(mantissa,divisor))); }

  // comparison of two Fixed-precision numbers
  bool operator> (Fixed<RADIX,EXPONENT> const other) const { return (mantissa >  other.mantissa); }
  bool operator< (Fixed<RADIX,EXPONENT> const other) const { return (mantissa <  other.mantissa); }
  bool operator>=(Fixed<RADIX,EXPONENT> const other) const { return (mantissa >= other.mantissa); }
  bool operator<=(Fixed<RADIX,EXPONENT> const other) const { return (mantissa <= other.mantissa); }
  bool operator==(Fixed<RADIX,EXPONENT> const other) const { return (mantissa == other.mantissa); }
  bool operator!=(Fixed<RADIX,EXPONENT> const other) const { return (mantissa != other.mantissa); }

private:
  
  // private constructor method for a Fixed-precision number given its mantissa
  Fixed(Integer n) : mantissa(n) { }
  
}; // Fixed

#endif // FIXED_H
