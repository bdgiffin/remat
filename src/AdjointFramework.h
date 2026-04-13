#ifndef ADJOINT_FRAMEWORK_H
#define ADJOINT_FRAMEWORK_H

#include "types.h"

enum class AdjointSupportLevel : int {
  Unsupported = 0,
  ScaffoldOnly = 1,
  Supported = 2
};

enum class ObjectivePolicyType : int {
  None = 0,
  TrussStressSquaredOverE = 1
};

struct ScalarStepInput {
  Real strain_n = 0.0;
  Real strain_np1 = 0.0;
  Real dt = 0.0;
};

struct ScalarStepRecord {
  Real dt = 0.0;
  Real strain_n = 0.0;
  Real strain_np1 = 0.0;
};

struct ScalarLocalSeed {
  Real bar_sigma_n = 0.0;
  Real bar_sigma_np1 = 0.0;
  Real direct_dE = 0.0;
};

struct ScalarGradientState {
  Real df_dtau = 0.0;
  Real df_dE = 0.0;
};

#endif // ADJOINT_FRAMEWORK_H
