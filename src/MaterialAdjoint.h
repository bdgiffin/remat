#ifndef MATERIAL_ADJOINT_H
#define MATERIAL_ADJOINT_H

#include <type_traits>
#include <utility>
#include "types.h"

namespace material_adjoint_detail {

template<typename M, typename = void>
struct has_num_params : std::false_type { };

template<typename M>
struct has_num_params<M, std::void_t<decltype(std::declval<const M&>().adjoint_num_params())> > : std::true_type { };

template<typename M, typename = void>
struct has_param_name : std::false_type { };

template<typename M>
struct has_param_name<M, std::void_t<decltype(std::declval<const M&>().adjoint_param_name(0))> > : std::true_type { };

template<typename M, typename = void>
struct has_scalar_stress_seed : std::false_type { };

template<typename M>
struct has_scalar_stress_seed<M, std::void_t<decltype(std::declval<M&>().adjoint_add_stress_seed(std::declval<Real*>(),Real(0.0)))> > : std::true_type { };

template<typename M, typename = void>
struct has_tensor_stress_seed : std::false_type { };

template<typename M>
struct has_tensor_stress_seed<M, std::void_t<decltype(std::declval<M&>().adjoint_add_stress_seed(std::declval<Real*>(),Real(0.0),Real(0.0),Real(0.0)))> > : std::true_type { };

template<typename M, typename = void>
struct has_objective_seed : std::false_type { };

template<typename M>
struct has_objective_seed<M, std::void_t<decltype(std::declval<M&>().adjoint_objective_seed(std::declval<Real*>()))> > : std::true_type { };

template<typename M, typename = void>
struct has_get_param_gradient : std::false_type { };

template<typename M>
struct has_get_param_gradient<M, std::void_t<decltype(std::declval<const M&>().adjoint_get_param_gradient(std::declval<const Real*>(),0))> > : std::true_type { };

template<typename M, typename = void>
struct has_clear_step_seed : std::false_type { };

template<typename M>
struct has_clear_step_seed<M, std::void_t<decltype(std::declval<M&>().adjoint_clear_step_seed(std::declval<Real*>()))> > : std::true_type { };

} // namespace material_adjoint_detail

template<typename M>
inline int material_adjoint_num_params(const M& model) {
  if constexpr (material_adjoint_detail::has_num_params<M>::value) {
    return int(model.adjoint_num_params());
  }
  return 0;
}

template<typename M>
inline const char* material_adjoint_param_name(const M& model, int i) {
  if constexpr (material_adjoint_detail::has_param_name<M>::value) {
    return model.adjoint_param_name(i);
  }
  (void)model;
  (void)i;
  return "";
}

template<typename M>
inline void material_adjoint_add_stress_seed(M& model, Real* state, Real seed) {
  if constexpr (material_adjoint_detail::has_scalar_stress_seed<M>::value) {
    model.adjoint_add_stress_seed(state,seed);
  } else {
    (void)model;
    (void)state;
    (void)seed;
  }
}

template<typename M>
inline void material_adjoint_add_stress_seed(M& model, Real* state, Real seed_xx, Real seed_yy, Real seed_xy) {
  if constexpr (material_adjoint_detail::has_tensor_stress_seed<M>::value) {
    model.adjoint_add_stress_seed(state,seed_xx,seed_yy,seed_xy);
  } else {
    (void)model;
    (void)state;
    (void)seed_xx;
    (void)seed_yy;
    (void)seed_xy;
  }
}

template<typename M>
inline void material_adjoint_objective_seed(M& model, Real* state) {
  if constexpr (material_adjoint_detail::has_objective_seed<M>::value) {
    model.adjoint_objective_seed(state);
  } else {
    (void)model;
    (void)state;
  }
}

template<typename M>
inline Real material_adjoint_get_param_gradient(const M& model, const Real* state, int i) {
  if constexpr (material_adjoint_detail::has_get_param_gradient<M>::value) {
    return model.adjoint_get_param_gradient(state,i);
  }
  (void)model;
  (void)state;
  (void)i;
  return Real(0.0);
}

template<typename M>
inline void material_adjoint_clear_step_seed(M& model, Real* state) {
  if constexpr (material_adjoint_detail::has_clear_step_seed<M>::value) {
    model.adjoint_clear_step_seed(state);
  } else {
    (void)model;
    (void)state;
  }
}

#endif // MATERIAL_ADJOINT_H
