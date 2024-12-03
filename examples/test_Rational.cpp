#include "Rational.h"
#include <iostream>

int main(void) {

  // Test creation of Rational numbers
  Rational a = Rational(5,6);
  Rational b = Rational(8,12);
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;

  // Test addition of two rational numbers
  auto c = a + b;
  std::cout << "a + b = " << c << std::endl;

  // Test subtraction of two rational numbers
  auto d = a - b;
  std::cout << "a - b = " << d << std::endl;
  
  return 0;
}
