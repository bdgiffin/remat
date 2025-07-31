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

  // NOTE: the sum and difference operations introduce the possibility for the sign of the first and second values to become different (try to eliminate this entirely!)
  
  // basic arithmetic operations between two Dual numbers of the same templated type
  Dual<T> operator+(Dual<T> const summand)    const {
    //Integer res;
    //if (__builtin_add_overflow(first.mantissa,summand.first.mantissa,&res))   std::cerr << "OVERFLOW: (+) " << first.mantissa  << " + " << summand.first.mantissa  << std::endl;
    //if (__builtin_add_overflow(second.mantissa,summand.second.mantissa,&res)) std::cerr << "OVERFLOW: (+) " << second.mantissa << " + " << summand.second.mantissa << std::endl;
    Dual<T> sum(first + summand.first,second + summand.second);
    // only do something if the sign of the first value changed (the first value must be non-zero
    //sum.second = copysign(sum.second,sum.first);
    return sum;
  }
  Dual<T> operator-(Dual<T> const subtrahend) const {
    //Integer res;
    //if (__builtin_sub_overflow(first.mantissa,subtrahend.first.mantissa,&res))   std::cerr << "OVERFLOW: (-) " << first.mantissa  << " - " << subtrahend.first.mantissa  << std::endl;
    //if (__builtin_sub_overflow(second.mantissa,subtrahend.second.mantissa,&res)) std::cerr << "OVERFLOW: (-) " << second.mantissa << " - " << subtrahend.second.mantissa << std::endl;
    Dual<T> difference(first - subtrahend.first,second - subtrahend.second);
    // only do something if the sign of the first value changed (the first value must be non-zero
    //difference.second = copysign(difference.second,difference.first);
    return difference;
  }

  // basic arithmetic operations between a Dual number and an Integer
  Dual<T> operator*(Integer const multiplier) const {
    LongInteger res;
    //if (__builtin_mul_overflow(first.mantissa,multiplier,&res)) std::cerr << "OVERFLOW: (*) " << first.mantissa << " * " << multiplier << std::endl;
    if ((first < zeroT) != (second < zeroT)) {
      return Dual<T>(first * multiplier - second % multiplier, second / multiplier);
    } else {
      return Dual<T>(first * multiplier + second % multiplier, second / multiplier);
    }
  }
  Dual<T> operator/(Integer const divisor)    const {
    LongInteger res;
    //if (__builtin_mul_overflow(second.mantissa,divisor,&res)) std::cerr << "OVERFLOW: (/) " << second.mantissa << " * " << divisor << std::endl;
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
    //LongInteger res;

    // Determine the sign of each mantissa
    bool first_signbit  = (long_first  < 0);
    bool second_signbit = (long_second < 0);

    // Conditionally shift the ranges for negative values
    if (first_signbit)  { long_first  += 1; }
    if (second_signbit) { long_second += 1; }
    
    // 1) apply multiplier's numerator
    //if (__builtin_mul_overflow(long_first,multiplier.numerator,&res)) std::cerr << "OVERFLOW: (*) " << long_first << " * " << multiplier.numerator << std::endl;
    //if ((first < zeroT) != (second < zeroT)) { // FIXME: there seems to be a corner case issue with the sign-dependent multiplication/division operation...?
    //  long_first  = long_first * multiplier.numerator - long_second % multiplier.numerator;
    //  long_second = long_second / multiplier.numerator;
    //} else {
    //  long_first  = long_first * multiplier.numerator + long_second % multiplier.numerator;
    //  long_second = long_second / multiplier.numerator;
    //}
    //if (first < zeroT) {
    if (first_signbit) {
      long_first  = long_first * multiplier.numerator - std::abs(long_second) % multiplier.numerator;
    } else {
      long_first  = long_first * multiplier.numerator + std::abs(long_second) % multiplier.numerator;
    }
    //long_first  = long_first  * multiplier.numerator + long_second % multiplier.numerator; // assuming first and second have same sign
    long_second = long_second / multiplier.numerator;
    
    //// 2) apply multiplier's denominator
    ////if (__builtin_mul_overflow(long_second,multiplier.denominator,&res)) std::cerr << "OVERFLOW: (*) " << long_second << " * " << multiplier.denominator << std::endl;
    //if ((first < zeroT) != (second < zeroT)) {
    //  long_second = long_second * multiplier.denominator - long_first % multiplier.denominator;
    //  long_first  = long_first / multiplier.denominator;
    //} else {
    //  long_second = long_second * multiplier.denominator + long_first % multiplier.denominator;
    //  long_first  = long_first / multiplier.denominator;
    //}
    //if (second < zeroT) {
    if (second_signbit) {
      long_second = long_second * multiplier.denominator - std::abs(long_first) % multiplier.denominator;
    } else {
      long_second = long_second * multiplier.denominator + std::abs(long_first) % multiplier.denominator;
    }
    //long_second = long_second * multiplier.denominator + long_first % multiplier.denominator; // assuming first and second have same sign
    long_first  = long_first  / multiplier.denominator;

    // Conditionally shift the ranges back for negative values
    if (first_signbit)  { long_first  -= 1; }
    if (second_signbit) { long_second -= 1; }

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
    //LongInteger res;

    // Determine the sign of each mantissa
    bool first_signbit  = (long_first  < 0);
    bool second_signbit = (long_second < 0);

    // Conditionally shift the ranges for negative values
    if (first_signbit)  { long_first  += 1; }
    if (second_signbit) { long_second += 1; }
    
    // first: apply divisor's denominator
    //if (__builtin_mul_overflow(long_first,divisor.denominator,&res)) std::cerr << "OVERFLOW: (*) " << long_first << " * " << divisor.denominator << std::endl;
    //if ((first < zeroT) != (second < zeroT)) {
    //  long_first  = long_first * divisor.denominator - long_second % divisor.denominator;
    //  long_second = long_second / divisor.denominator;
    //} else {
    //  long_first  = long_first * divisor.denominator + long_second % divisor.denominator;
    //  long_second = long_second / divisor.denominator;
    //}
    //if (first < zeroT) {
    if (first_signbit) {
      long_first  = long_first * divisor.denominator - std::abs(long_second) % divisor.denominator;
    } else {
      long_first  = long_first * divisor.denominator + std::abs(long_second) % divisor.denominator;
    }
    //long_first  = long_first  * divisor.denominator + long_second % divisor.denominator; // assuming first and second have same sign
    long_second = long_second / divisor.denominator;
    
    // second: apply divisor's numerator
    //if (__builtin_mul_overflow(long_second,divisor.numerator,&res)) std::cerr << "OVERFLOW: (*) " << long_second << " * " << divisor.numerator << std::endl;
    //if ((first < zeroT) != (second < zeroT)) {
    //  long_second = long_second * divisor.numerator - long_first % divisor.numerator;
    //  long_first  = long_first / divisor.numerator;
    //} else {
    //  long_second = long_second * divisor.numerator + long_first % divisor.numerator;
    //  long_first  = long_first / divisor.numerator;
    //}
    //if (second < zeroT) {
    if (second_signbit) {
      long_second = long_second * divisor.numerator - std::abs(long_first) % divisor.numerator;
    } else {
      long_second = long_second * divisor.numerator + std::abs(long_first) % divisor.numerator;
    }
    //long_second = long_second * divisor.numerator + long_first % divisor.numerator; // assuming first and second have same sign
    long_first  = long_first  / divisor.numerator;
    
    // Conditionally shift the ranges back for negative values
    if (first_signbit)  { long_first  -= 1; }
    if (second_signbit) { long_second -= 1; }

    // downcast and return result
    T new_first, new_second;
    new_first.mantissa  = long_first;
    new_second.mantissa = long_second;
    return Dual<T>(new_first,new_second);
  }

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
