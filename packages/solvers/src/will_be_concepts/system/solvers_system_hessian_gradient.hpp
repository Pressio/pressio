

#ifndef SOLVERS_SYSTEM_MEETS_HESSIAN_GRADIENT_API_HPP_
#define SOLVERS_SYSTEM_MEETS_HESSIAN_GRADIENT_API_HPP_

namespace pressio{ namespace solvers{ namespace concepts {

template<typename T, typename enable = void>
struct system_hessian_gradient : std::false_type{};

template<typename T>
struct system_hessian_gradient
<T,
 ::pressio::mpl::enable_if_t<
   ::pressio::solvers::predicates::has_scalar_typedef<T>::value   and
   ::pressio::solvers::predicates::has_state_typedef<T>::value    and
   ::pressio::solvers::predicates::has_hessian_typedef<T>::value and
   ::pressio::solvers::predicates::has_gradient_typedef<T>::value and

   ::pressio::solvers::predicates::has_const_create_hessian_method_return_result<
      T, typename T::hessian_type>::value and

   ::pressio::solvers::predicates::has_const_create_gradient_method_return_result<
      T, typename T::gradient_type>::value and

   ::pressio::solvers::predicates::has_const_hessian_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::hessian_type>::value and

   ::pressio::solvers::predicates::has_const_gradient_method_accept_state_result_norm_return_void<
      T, typename T::state_type, typename T::gradient_type, typename T::scalar_type>::value 
   >
 > : std::true_type{};

}}} // namespace pressio::solvers::concepts
#endif
