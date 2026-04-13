#ifndef CONSTITUTIVE_ADJOINT_H
#define CONSTITUTIVE_ADJOINT_H

#include <string>
#include <vector>

template<class RealT>
struct MaterialPointKinematics {
  // Generic scalar kinematic measure; interpretation is model-specific.
  RealT primary_measure = RealT(0.0);
  RealT dt = RealT(0.0);
};

template<class RealT>
struct MaterialForwardOutput {
  RealT stress = RealT(0.0);
  RealT psi = RealT(0.0);
};

template<class RealT>
struct MaterialReverseSeed {
  RealT dL_dstress = RealT(0.0);
};

template<class RealT>
struct MaterialReverseOutput {
  RealT dL_dmeasure = RealT(0.0);
  RealT stress = RealT(0.0);
  RealT psi = RealT(0.0);
  std::vector<RealT> dL_dparams;
};

template<class RealT>
class ParameterAccumulator {
public:
  void clear() {
    m_param_names.clear();
    m_param_grads.clear();
  }

  int register_param(const std::string& name) {
    for (int i=0; i<int(m_param_names.size()); ++i) {
      if (m_param_names[i] == name) { return i; }
    }
    m_param_names.push_back(name);
    m_param_grads.push_back(RealT(0.0));
    return int(m_param_names.size()) - 1;
  }

  void zero() {
    for (int i=0; i<int(m_param_grads.size()); ++i) {
      m_param_grads[i] = RealT(0.0);
    }
  }

  int size() const { return int(m_param_names.size()); }

  const std::string& name(int param_id) const { return m_param_names[param_id]; }

  RealT get(int param_id) const { return m_param_grads[param_id]; }

  void set(int param_id, RealT value) {
    if (param_id < 0 || param_id >= int(m_param_grads.size())) { return; }
    m_param_grads[param_id] = value;
  }

  void add(int param_id, RealT value) {
    if (param_id < 0 || param_id >= int(m_param_grads.size())) { return; }
    m_param_grads[param_id] += value;
  }

private:
  std::vector<std::string> m_param_names;
  std::vector<RealT> m_param_grads;
};

#endif // CONSTITUTIVE_ADJOINT_H
