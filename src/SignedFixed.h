#ifndef FIXED_H
#define FIXED_H

#include "types.h"
#include <iostream>
#include <math.h>
#include <cmath>

template<int RADIX = 10, int EXPONENT = 0>
struct Fixed {
  inline static Real     decimal = std::pow(RADIX,+EXPONENT);
  inline static Real inv_decimal = std::pow(RADIX,-EXPONENT);
  Integer mantissa;
  bool    signbit;  // negative if true (positive if false)

  // constructor method for a Fixed-precision number
  Fixed(double d = 0.0) : mantissa(std::abs(d*inv_decimal)), signbit(std::signbit(d)) { }
  
  // private constructor method for a Fixed-precision number given its mantissa and signbit
  Fixed(Integer n, bool pm) : mantissa(std::abs(n)), signbit(pm) { }

  // Convert to a signed integer
  Integer to_Integer(void) const { return signbit ? -mantissa : mantissa; }

  // implicit type-cast to a floating point number
  operator Real() const { return to_Integer()*decimal; }

  // WARNING: the summation operation is no longer commutative! (but it is reversible?)
  // basic arithmetic operations between two Fixed-precision numbers with the same radix and exponent
  Fixed<RADIX,EXPONENT> operator+(Fixed<RADIX,EXPONENT> const summand)    const {
    Integer sum = to_Integer() + summand.to_Integer();
    if (sum == 0) {
      return Fixed<RADIX,EXPONENT>(0,signbit);
    } else {
      return Fixed<RADIX,EXPONENT>(sum);
    }
  }
  Fixed<RADIX,EXPONENT> operator-(Fixed<RADIX,EXPONENT> const subtrahend) const {
    Integer difference = to_Integer() - subtrahend.to_Integer();
    if (difference == 0) {
      return Fixed<RADIX,EXPONENT>(0,signbit);
    } else {
      return Fixed<RADIX,EXPONENT>(difference);
    }
  }

  // basic arithmetic operations between a Fixed-precision number and an Integer
  Fixed<RADIX,EXPONENT> operator*(Integer const multiplier) const { return Fixed<RADIX,EXPONENT>(to_Integer()*multiplier); }
  Fixed<RADIX,EXPONENT> operator/(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(to_Integer()/divisor);    }
  Fixed<RADIX,EXPONENT> operator%(Integer const divisor)    const { return Fixed<RADIX,EXPONENT>(to_Integer()%divisor);    }

  // basic arithmetic operations between a Fixed-precision number and a Real (floating point) number
  Fixed<RADIX,EXPONENT> operator*(Real const multiplier) const { return Fixed<RADIX,EXPONENT>(to_Integer()*multiplier); }
  Fixed<RADIX,EXPONENT> operator/(Real const divisor)    const { return Fixed<RADIX,EXPONENT>(to_Integer()/divisor);    }
  Fixed<RADIX,EXPONENT> operator%(Real const divisor)    const { return Fixed<RADIX,EXPONENT>(std::remainder(to_Integer(),divisor)); }

  // comparison of two Fixed-precision numbers
  bool operator> (Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() >  other.to_Integer()); }
  bool operator< (Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() <  other.to_Integer()); }
  bool operator>=(Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() >= other.to_Integer()); }
  bool operator<=(Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() <= other.to_Integer()); }
  bool operator==(Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() == other.to_Integer()); }
  bool operator!=(Fixed<RADIX,EXPONENT> const other) const { return (to_Integer() != other.to_Integer()); }

private:
  
  // private constructor method for a Fixed-precision number given a signed Integer
  Fixed(Integer n) : mantissa(std::abs(n)), signbit(std::signbit(n)) { }
  
}; // Fixed

// Declare standard Fixed-precision numbers
typedef Fixed<10,-6> Fixed_V;
typedef Fixed<10,-6> Fixed_U;
typedef Fixed<10,-6> Fixed_E;

// Declare conversion to/from a Real value
void load_from_Real(Real value, Fixed_E &load_value) {
  //std::memcpy(&(load_value.mantissa),&value,sizeof(Integer));
  load_value.mantissa = std::abs(Integer(value));
  load_value.signbit  = std::signbit(value);
}
void save_as_Real(Fixed_E value, Real& save_value) {
  //std::memcpy(&save_value,&(value.mantissa),sizeof(Integer));
  save_value = value.to_Integer();
}


#endif // FIXED_H
