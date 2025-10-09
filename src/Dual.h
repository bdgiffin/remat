#ifndef DUAL_H
#define DUAL_H

#include "types.h"
#include "arithmetic.h"
#include "Rational.h"
#include <iostream>
#include <math.h>
#include <tuple>

template<class T>
struct Dual {
  T first, second;

  // constructor for a Dual number
  Dual(T a, T b = zeroT) : first(a), second(b) { }

  // implicit type-cast to a primal number
  operator T() const { return first; }
  
  // basic arithmetic operations between two Dual numbers of the same templated type
  Dual<T> operator+(Dual<T> const summand)    const {
    //Integer res;
    //if (__builtin_add_overflow(first.mantissa,summand.first.mantissa,&res))   std::cerr << "OVERFLOW: (+) " << first.mantissa  << " + " << summand.first.mantissa  << std::endl;
    //if (__builtin_add_overflow(second.mantissa,summand.second.mantissa,&res)) std::cerr << "OVERFLOW: (+) " << second.mantissa << " + " << summand.second.mantissa << std::endl;
    Dual<T> sum(first + summand.first,second + summand.second);
    return sum;
  }
  Dual<T> operator-(Dual<T> const subtrahend) const {
    //Integer res;
    //if (__builtin_sub_overflow(first.mantissa,subtrahend.first.mantissa,&res))   std::cerr << "OVERFLOW: (-) " << first.mantissa  << " - " << subtrahend.first.mantissa  << std::endl;
    //if (__builtin_sub_overflow(second.mantissa,subtrahend.second.mantissa,&res)) std::cerr << "OVERFLOW: (-) " << second.mantissa << " - " << subtrahend.second.mantissa << std::endl;
    Dual<T> difference(first - subtrahend.first,second - subtrahend.second);
    return difference;
  }

  // basic arithmetic operations between a Dual number and an Integer
  Dual<T> operator*(Integer const multiplier) const { return Dual<T>(first * multiplier + second % multiplier, second / multiplier); }
  Dual<T> operator/(Integer const divisor)    const { return Dual<T>(first / divisor, second * divisor + first % divisor); }

  // basic arithmetic operations between a Dual number and a Rational number (may overflow!)
  Dual<T> operator*(Rational multiplier) const {
    T new_first, new_second;
    auto [a, b] = squeezeE(first.mantissa, second.mantissa, multiplier.numerator, multiplier.denominator);
    new_first.mantissa  = a;
    new_second.mantissa = b;
    return Dual<T>(new_first,new_second);
  } // operator*(Rational)
  Dual<T> operator/(Rational divisor)    const {
    T new_first, new_second;
    auto [a, b] = squeezeE(first.mantissa, second.mantissa, divisor.denominator, divisor.numerator);
    new_first.mantissa  = a;
    new_second.mantissa = b;
    return Dual<T>(new_first,new_second);
  } // operator/(Rational)

  // basic arithmetic operations between a Dual number and a Real (floating point) number
  Dual<T> operator*(Real const multiplier) const { return Dual<T>(first*multiplier,second/multiplier); }
  Dual<T> operator/(Real const divisor)    const { return Dual<T>(first/divisor,second*divisor); }

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
