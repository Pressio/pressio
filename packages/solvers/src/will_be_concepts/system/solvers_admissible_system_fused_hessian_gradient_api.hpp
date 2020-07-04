

#ifndef SOLVERS_SYSTEM_MEETS_FUSED_HESSIAN_GRADIENT_API_HPP_
#define SOLVERS_SYSTEM_MEETS_FUSED_HESSIAN_GRADIENT_API_HPP_

namespace pressio{ namespace solvers{ namespace meta {

template<typename T, typename enable = void>
struct system_meets_fused_hessian_gradient_api : std::false_type{};

template<typename T>
struct system_meets_fused_hessian_gradient_api
<T,
 ::pressio::mpl::enable_if_t<
   ::pressio::solvers::meta::has_scalar_typedef<T>::value   and
   ::pressio::solvers::meta::has_state_typedef<T>::value    and
   ::pressio::solvers::meta::has_hessian_typedef<T>::value and
   ::pressio::solvers::meta::has_gradient_typedef<T>::value and

   ::pressio::solvers::meta::has_const_create_hessian_method_return_result<
      T, typename T::hessian_type>::value and

   ::pressio::solvers::meta::has_const_create_gradient_method_return_result<
      T, typename T::gradient_type>::value and

   ::pressio::solvers::meta::has_const_hessianandgradient_method_accept_state_result_norm_return_void<
      T, typename T::state_type, 
      typename T::hessian_type, typename T::gradient_type, 
      typename T::scalar_type>::value and

   ::pressio::solvers::meta::has_const_residualnorm_method_accept_state_norm_return_void<
      T, typename T::state_type, typename T::scalar_type>::value 
   >
 > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
