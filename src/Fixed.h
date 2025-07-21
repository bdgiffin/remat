#ifndef FIXED_H
#define FIXED_H

#include "types.h"
#include <iostream>
#include <math.h>

template<int RADIX = 10, int EXPONENT = 0>
struct Fixed {
  inline static Real     decimal = std::pow(RADIX,+EXPONENT);
  inline static Real inv_decimal = std::pow(RADIX,-EXPONENT);
  Integer mantissa;

  // constructor method for a Fixed-precision number
  Fixed(double d = 0.0) : mantissa(d*inv_decimal) { }

  // implicit type-cast to a floating point number
  operator Real() const { return mantissa*decimal; }

  // basic arithmetic operations between two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator+(Fixed<RADIX,EXPONENT> const summand)    const { return Fixed<RADIX,EXPONENT>(mantissa + summand.mantissa);    }
  Fixed<RADIX,EXPONENT> operator-(Fixed<RADIX,EXPONENT> const subtrahend) const { return Fixed<RADIX,EXPONENT>(mantissa - subtrahend.mantissa); }

  // basic arithmetic operations between a Fixed-precision number and an Integer
  Fixed<RADIX,EXPONENT> operator*(Integer const multiplier) const { return Fixed<RADIX,EXPONENT>(mantissa*multiplier); }
  Fixed<RADIX,EXPONENT> operator/(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(mantissa/divisor);    }
  Fixed<RADIX,EXPONENT> operator%(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(mantissa%divisor);    }

  // basic arithmetic operations between a Fixed-precision number and a Real (floating point) number
  Fixed<RADIX,EXPONENT> operator*(Real const multiplier) const { return Fixed<RADIX,EXPONENT>(Integer(mantissa*multiplier)); }
  Fixed<RADIX,EXPONENT> operator/(Real const divisor)    const { return Fixed<RADIX,EXPONENT>(Integer(mantissa/divisor));    }
  Fixed<RADIX,EXPONENT> operator%(Real const divisor)    const { return Fixed<RADIX,EXPONENT>(Integer(std::remainder(mantissa,divisor))); }

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

// Declare standard Fixed-precision numbers
typedef Fixed<10,-6> Fixed_V;
typedef Fixed<10,-6> Fixed_U;
typedef Fixed<10,-6> Fixed_E;

// Declare conversion to/from a Real value
void load_from_Real(Real value, Fixed_E &load_value) { 
  std::memcpy(&(load_value.mantissa),&value,sizeof(Integer));
}
void save_as_Real(Fixed_E value, Real& save_value) {
  std::memcpy(&save_value,&(value.mantissa),sizeof(Integer));
}


#endif // FIXED_H
