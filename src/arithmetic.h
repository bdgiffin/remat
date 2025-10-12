#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include "types.h"
#include <tuple>

// Implementations for Euclidean and floored division adapted from: 
// "Division and Modulus for Computer Scientists" - Daan Leijen (2001)
// https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/divmodnote-letter.pdf

// Truncated division
template<class T>
inline std::tuple<T, T> divmodT(T dividend, T divisor) {
  T quotient  = dividend/divisor;
  T remainder = dividend%divisor; // dividend - quotient*divisor
  return {quotient, remainder};
} // divmodT

// Euclidean division
template<class T>
inline std::tuple<T, T> divmodE(T dividend, T divisor) {
  auto [quotient, remainder] = divmodT(dividend, divisor);
  if (remainder < 0) {
    if (divisor > 0) {
      quotient  -= 1;
      remainder += divisor;
    } else {
      quotient  += 1;
      remainder -= divisor;
    }
  }
  return {quotient, remainder};
} // divmodE
  
// Floored division
template<class T>
inline std::tuple<T, T> divmodF(T dividend, T divisor) {
  auto [quotient, remainder] = divmodT(dividend, divisor);
  if ((remainder > 0 && divisor < 0) || (remainder < 0 && divisor > 0)) {
    quotient  -= 1;
    remainder += divisor;
  }
  return {quotient, remainder};
} // divmodF

// Discrete squeeze mapping involving a rational squeeze parameter (p/q) based on truncated division
inline std::tuple<Integer, Integer> squeezeT(LongInteger first, LongInteger second, LongInteger p, LongInteger q) {
  // ASSUMES BOTH P AND Q ARE NON-NEGATIVE!
  
  // Determine the sign of the first and second arguments
  bool first_signbit  = (first  < 0);
  bool second_signbit = (second < 0);

  // Conditionally shift the ranges for negative values
  if (first_signbit)  { first  += 1; }
  if (second_signbit) { second += 1; }
    
  // step 1) squeeze mapping using the squeeze parameter p
  //if (__builtin_mul_overflow(first,p,&res)) std::cerr << "OVERFLOW: (*) " << first << " * " << p << std::endl;
  if (first_signbit) {
    first = first * p - std::abs(second) % p;
  } else {
    first = first * p + std::abs(second) % p;
  }
  second = second / p;
    
  // step 2) inverse squeeze mapping using the squeeze parameter q
  //if (__builtin_mul_overflow(second,q,&res)) std::cerr << "OVERFLOW: (*) " << second << " * " << q << std::endl;
  if (second_signbit) {
    second = second * q - std::abs(first) % q;
  } else {
    second = second * q + std::abs(first) % q;
  }
  first = first / q;

  // Conditionally shift the ranges back for negative values
  if (first_signbit)  { first  -= 1; }
  if (second_signbit) { second -= 1; }

  // Downcast and return result
  return {Integer(first), Integer(second)};
} // squeezeT

// Discrete squeeze mapping involving a rational squeeze parameter (p/q) based on Euclidean division
inline std::tuple<Integer, Integer> squeezeE(LongInteger first, LongInteger second, LongInteger p, LongInteger q) {
  // ASSUMES BOTH P AND Q ARE NON-NEGATIVE!
    
  // step 1) squeeze mapping using the squeeze parameter p
  //if (__builtin_mul_overflow(first,p,&res)) std::cerr << "OVERFLOW: (*) " << first << " * " << p << std::endl;
  auto [quotient1, remainder1] = divmodE(second, p);
  first  = first * p + remainder1;
  second = quotient1;
    
  // step 2) inverse squeeze mapping using the squeeze parameter q
  //if (__builtin_mul_overflow(second,q,&res)) std::cerr << "OVERFLOW: (*) " << second << " * " << q << std::endl;
  auto [quotient2, remainder2] = divmodE(first, q);
  second = second * q + remainder2;
  first  = quotient2;

  // Downcast and return result
  return {Integer(first), Integer(second)};
} // squeezeE

// Discrete squeeze mapping involving a rational squeeze parameter (p/q) based on floored division
inline std::tuple<Integer, Integer> squeezeF(LongInteger first, LongInteger second, LongInteger p, LongInteger q) {
  // ASSUMES BOTH P AND Q ARE NON-NEGATIVE!
    
  // step 1) squeeze mapping using the squeeze parameter p
  //if (__builtin_mul_overflow(first,p,&res)) std::cerr << "OVERFLOW: (*) " << first << " * " << p << std::endl;
  auto [quotient1, remainder1] = divmodF(second, p);
  first  = first * p + remainder1;
  second = quotient1;
    
  // step 2) inverse squeeze mapping using the squeeze parameter q
  //if (__builtin_mul_overflow(second,q,&res)) std::cerr << "OVERFLOW: (*) " << second << " * " << q << std::endl;
  auto [quotient2, remainder2] = divmodF(first, q);
  second = second * q + remainder2;
  first  = quotient2;

  // Downcast and return result
  return {Integer(first), Integer(second)};
} // squeezeF

#endif // ARITHMETIC_H
