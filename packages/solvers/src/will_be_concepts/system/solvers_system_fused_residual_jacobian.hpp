
#ifndef SOLVERS_SYSTEM_FUSED_MEETS_RESIDUAL_JACOBIAN_API_HPP_
#define SOLVERS_SYSTEM_FUSED_MEETS_RESIDUAL_JACOBIAN_API_HPP_

namespace pressio{ namespace solvers{ namespace concepts {

template<typename T, typename enable = void>
struct system_fused_residual_jacobian : std::false_type{};

template<typename T>
struct system_fused_residual_jacobian
<T,
 ::pressio::mpl::enable_if_t<
   ::pressio::solvers::predicates::has_scalar_typedef<T>::value   and
   ::pressio::solvers::predicates::has_state_typedef<T>::value    and
   ::pressio::solvers::predicates::has_residual_typedef<T>::value and
   ::pressio::solvers::predicates::has_jacobian_typedef<T>::value and

   ::pressio::solvers::predicates::has_const_create_residual_method_return_result<
      T, typename T::residual_type>::value and

   ::pressio::solvers::predicates::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type>::value and

   ::pressio::solvers::predicates::has_const_residualandjacobian_method_accept_state_result_norm_return_void<
      T, typename T::state_type, typename T::residual_type, 
      typename T::jacobian_type, typename T::scalar_type>::value and

   ::pressio::solvers::predicates::has_const_residualnorm_method_accept_state_norm_return_void<
      T, typename T::state_type, typename T::scalar_type>::value    
  >
 > : std::true_type{};

}}} // namespace pressio::solvers::concepts
#endif
