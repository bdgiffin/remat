#ifndef DUAL_H
#define DUAL_H

#include "types.h"
#include <iostream>

template<class T>
struct Dual {
  T first, second;

  // constructor for a Dual number
  Dual(T a, T b = 0) : first(a), second(b) { }

  // multiplication of a Dual number by an Integer
  Dual<T> operator*(Integer const multiplier) {
    return Dual<T>(first * multiplier + second % multiplier, second / multiplier);
  } // operator*

  // division of a Dual number by an Integer
  Dual<T> operator/(Integer const divisor) {
    return Dual<T>(first / divisor, second * divisor + first % divisor);
  } // operator/
  
}; // Dual

#endif // DUAL_H
