#include "Fixed.h"
#include <iostream>

int main(void) {

  // Test creation of Fixed-precision numbers
  typedef Fixed<10,-6> FixedE6;
  FixedE6 a = 5.0/6.0;
  FixedE6 b = 8.0/12.0;
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;
    
  // Test addition of two Fixed-precision numbers
  auto c = a + b;
  std::cout << "c = a + b = " << c << std::endl;

  // Test subtraction of two Fixed-precision numbers
  auto d = c - b;
  std::cout << "a = c - b = " << d << std::endl;
  
  return 0;
}
