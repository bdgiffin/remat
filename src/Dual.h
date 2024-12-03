#ifndef DUAL_H
#define DUAL_H

#include "types.h"
#include <iostream>

template<class T>
struct Dual {
  T first, second;

  // constructor for a Dual number
  Dual(T a, T b = 0) : first(a), second(b) { }
  
  // basic arithmetic operations between two Dual numbers of the same templated type
  Dual<T> operator+(Dual<T> const summand)    { return Dual<T>(first +    summand.first,second +    summand.second); }
  Dual<T> operator-(Dual<T> const subtrahend) { return Dual<T>(first - subtrahend.first,second - subtrahend.second); }

  // basic arithmetic operations between a Dual number and an Integer
  Dual<T> operator*(Integer const multiplier) { return Dual<T>(first * multiplier + second % multiplier, second / multiplier               ); }
  Dual<T> operator/(Integer const divisor)    { return Dual<T>(first / divisor                         , second * divisor + first % divisor); }

  // comparison of two Dual numbers
  bool operator==(Dual<T> const other) { return ((first == other.first) && (second == other.second)); }
  bool operator!=(Dual<T> const other) { return ((first != other.first) || (second != other.second)); }
  
}; // Dual

#endif // DUAL_H
