#ifndef DUAL_H
#define DUAL_H

#include "types.h"
#include "Rational.h"
#include <iostream>
#include <math.h>

template<class T>
struct Dual {
  T first, second;

  // constructor for a Dual number
  Dual(T a, T b = 0) : first(a), second(b) { }
  
  // basic arithmetic operations between two Dual numbers of the same templated type
  Dual<T> operator+(Dual<T> const summand)    const { return Dual<T>(first +    summand.first,second +    summand.second); }
  Dual<T> operator-(Dual<T> const subtrahend) const { return Dual<T>(first - subtrahend.first,second - subtrahend.second); }

  // basic arithmetic operations between a Dual number and an Integer
  Dual<T> operator*(Integer const multiplier) const {
    if ((first < zeroT) != (second < zeroT)) {
      return Dual<T>(first * multiplier - second % multiplier, second / multiplier);
    } else {
      return Dual<T>(first * multiplier + second % multiplier, second / multiplier);
    }
  }
  Dual<T> operator/(Integer const divisor)    const {
    if ((first < zeroT) != (second < zeroT)) {
      return Dual<T>(first / divisor, second * divisor - first % divisor);
    } else {
      return Dual<T>(first / divisor, second * divisor + first % divisor);
    }
  }

  // basic arithmetic operations between a Dual number and a Rational number (may overflow!)
  Dual<T> operator*(Rational const multiplier) const { return ((*this)*multiplier.numerator)/multiplier.denominator; }
  Dual<T> operator/(Rational const divisor)    const { return ((*this)*divisor.denominator)/divisor.numerator;       }

  // basic arithmetic operations between a Dual number and a double-precision floating point number
  Dual<T> operator*(double const multiplier) const { return (*this)*multiplier; }
  Dual<T> operator/(double const divisor)    const { return (*this)*divisor;    }

  // comparison of two Dual numbers
  bool operator==(Dual<T> const other) const { return ((first == other.first) && (second == other.second)); }
  bool operator!=(Dual<T> const other) const { return ((first != other.first) || (second != other.second)); }

private:
  inline static T zeroT = T(0.0);
}; // Dual

// overloaded output stream operator
template<class T>
std::ostream& operator<<(std::ostream& os, const Dual<T>& number) {
  os << '{' << number.first << ',' << number.second << '}';
  return os;
} // operator<<

#endif // DUAL_H
