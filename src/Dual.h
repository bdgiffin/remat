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
  Dual(T a, T b = zeroT) : first(a), second(b) { }

  // implicit type-cast to a primal number
  operator T() const { return first; }
  
  // basic arithmetic operations between two Dual numbers of the same templated type
  Dual<T> operator+(Dual<T> const summand)    const {
    //Integer res;
    //if (__builtin_add_overflow(first.mantissa,summand.first.mantissa,&res))   std::cerr << "OVERFLOW: (+) " << first.mantissa  << " + " << summand.first.mantissa  << std::endl;
    //if (__builtin_add_overflow(second.mantissa,summand.second.mantissa,&res)) std::cerr << "OVERFLOW: (+) " << second.mantissa << " + " << summand.second.mantissa << std::endl;
    return Dual<T>(first +    summand.first,second +    summand.second);
  }
  Dual<T> operator-(Dual<T> const subtrahend) const {
    //Integer res;
    //if (__builtin_sub_overflow(first.mantissa,subtrahend.first.mantissa,&res))   std::cerr << "OVERFLOW: (-) " << first.mantissa  << " - " << subtrahend.first.mantissa  << std::endl;
    //if (__builtin_sub_overflow(second.mantissa,subtrahend.second.mantissa,&res)) std::cerr << "OVERFLOW: (-) " << second.mantissa << " - " << subtrahend.second.mantissa << std::endl;
    return Dual<T>(first - subtrahend.first,second - subtrahend.second);
  }

  // basic arithmetic operations between a Dual number and an Integer
  Dual<T> operator*(Integer const multiplier) const {
    LongInteger res;
    if (__builtin_mul_overflow(first.mantissa,multiplier,&res)) std::cerr << "OVERFLOW: (*) " << first.mantissa << " * " << multiplier << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      return Dual<T>(first * multiplier - second % multiplier, second / multiplier);
    } else {
      return Dual<T>(first * multiplier + second % multiplier, second / multiplier);
    }
  }
  Dual<T> operator/(Integer const divisor)    const {
    LongInteger res;
    if (__builtin_mul_overflow(second.mantissa,divisor,&res)) std::cerr << "OVERFLOW: (/) " << second.mantissa << " * " << divisor << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      return Dual<T>(first / divisor, second * divisor - first % divisor);
    } else {
      return Dual<T>(first / divisor, second * divisor + first % divisor);
    }
  }

  // basic arithmetic operations between a Dual number and a Rational number (may overflow!)
  //Dual<T> operator*(Rational const multiplier) const { return ((*this)*multiplier.numerator)/multiplier.denominator; }
  //Dual<T> operator/(Rational const divisor)    const { return ((*this)*divisor.denominator)/divisor.numerator;       }

  Dual<T> operator*(Rational const multiplier) const {
    //return ((*this)*multiplier.numerator)/multiplier.denominator;

    // temporarily cast first and second mantissas to long integers
    LongInteger long_first  = first.mantissa;
    LongInteger long_second = second.mantissa;
    LongInteger res;
    
    // 1) apply multiplier's numerator
    if (__builtin_mul_overflow(long_first,multiplier.numerator,&res)) std::cerr << "OVERFLOW: (*) " << long_first << " * " << multiplier.numerator << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      long_first  = long_first * multiplier.numerator - long_second % multiplier.numerator;
      long_second = long_second / multiplier.numerator;
    } else {
      long_first  = long_first * multiplier.numerator + long_second % multiplier.numerator;
      long_second = long_second / multiplier.numerator;
    }
    
    // 2) apply multiplier's denominator
    if (__builtin_mul_overflow(long_second,multiplier.denominator,&res)) std::cerr << "OVERFLOW: (*) " << long_second << " * " << multiplier.denominator << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      long_second = long_second * multiplier.denominator - long_first % multiplier.denominator;
      long_first  = long_first / multiplier.denominator;
    } else {
      long_second = long_second * multiplier.denominator + long_first % multiplier.denominator;
      long_first  = long_first / multiplier.denominator;
    }

    // downcast and return result
    T new_first, new_second;
    new_first.mantissa  = long_first;
    new_second.mantissa = long_second;
    return Dual<T>(new_first,new_second);
  }

  Dual<T> operator/(Rational const divisor)    const {
    //return ((*this)*divisor.denominator)/divisor.numerator;
    
    // temporarily cast first and second mantissas to long integers
    LongInteger long_first  = first.mantissa;
    LongInteger long_second = second.mantissa;
    LongInteger res;
    
    // first: apply divisor's denominator
    if (__builtin_mul_overflow(long_first,divisor.denominator,&res)) std::cerr << "OVERFLOW: (*) " << long_first << " * " << divisor.denominator << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      long_first  = long_first * divisor.denominator - long_second % divisor.denominator;
      long_second = long_second / divisor.denominator;
    } else {
      long_first  = long_first * divisor.denominator + long_second % divisor.denominator;
      long_second = long_second / divisor.denominator;
    }
    
    // second: apply divisor's numerator
    if (__builtin_mul_overflow(long_second,divisor.numerator,&res)) std::cerr << "OVERFLOW: (*) " << long_second << " * " << divisor.numerator << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      long_second = long_second * divisor.numerator - long_first % divisor.numerator;
      long_first  = long_first / divisor.numerator;
    } else {
      long_second = long_second * divisor.numerator + long_first % divisor.numerator;
      long_first  = long_first / divisor.numerator;
    }

    // downcast and return result
    T new_first, new_second;
    new_first.mantissa  = long_first;
    new_second.mantissa = long_second;
    return Dual<T>(new_first,new_second);
  }

  // basic arithmetic operations between a Dual number and a double-precision floating point number
  Dual<T> operator*(double const multiplier) const { return Dual<T>(first*multiplier,second/multiplier); }
  Dual<T> operator/(double const divisor)    const { return Dual<T>(first/divisor,second*divisor); }

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
