#ifndef TYPES_H
#define TYPES_H

#include <cstdint>

// Define a Fixed-precision numeric type
typedef int32_t Integer;
typedef int64_t LongInteger;

// Define a standard Floating point ("Real") numeric type
typedef double Real;

// Declare conversion to/from a Real value
void load_from_Real(Real value, Real &load_value) { load_value = value; }
void save_as_Real(Real value, Real& save_value) { save_value = value; }

#endif // TYPES_H
