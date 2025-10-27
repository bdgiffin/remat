#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <cmath>
#include <limits>
#include <iostream>

// Define a Fixed-precision numeric type
typedef int32_t Integer;
typedef int64_t LongInteger;

// Define a standard Floating point ("Real") numeric type
typedef double Real;

// Declare conversion to/from a Real value
void load_from_Real(Real value, Real &load_value) { load_value = value; }
void save_as_Real(Real value, Real& save_value) { save_value = value; }

// Declare function to return the smallest representable value with the sign of the incoming argument
Real smallest_value(Real signed_value) { 
    std::cout << "USING TYPES_H smallest_value" << std::endl;
    return std::copysign(std::numeric_limits<Real>::min(),signed_value); 
}

// Return a value with the magnitude of the first argument and the sign of the second argument
Real copysign(Real magnitude, Real signed_value) { 
    std::cout << "USING TYPES_H copysign" << std::endl;
    return std::copysign(magnitude,signed_value); 
}

#endif // TYPES_H
