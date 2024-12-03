#include "Fixed.h"
#include "Dual.h"
#include <iostream>

int main(void) {

  // Test creation of Dual numbers
  typedef Fixed<10,-6> FixedE6;
  typedef Dual<FixedE6> DualE6;
  FixedE6 a = 5.0/6.0;
  FixedE6 b = 0.0;
  DualE6  d = DualE6(a,b);
  std::cout << "d.first  = " << d.first  << std::endl;
  std::cout << "d.second = " << d.second << std::endl;
    
  // Test multiplication of dual numbers
  auto q = d / 100;
  std::cout << "q.first  = " << q.first  << std::endl;
  std::cout << "q.second = " << q.second << std::endl;

  // Test subtraction of two Fixed-precision numbers
  auto m = q * 100;
  std::cout << "m.first  = " << m.first  << std::endl;
  std::cout << "m.second = " << m.second << std::endl;
  
  return 0;
}
