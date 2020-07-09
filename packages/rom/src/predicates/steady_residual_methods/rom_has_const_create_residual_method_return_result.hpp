
#ifndef rom_has_const_create_residual_method_return_result_HPP_
#define rom_has_const_create_residual_method_return_result_HPP_

namespace pressio{ namespace rom{ namespace predicates {
  
template <
  typename T,
  typename result_t,
  typename = void
  >
struct has_const_create_residual_method_return_result
  : std::false_type{};

template <
  typename T,
  typename result_t
  >
struct has_const_create_residual_method_return_result<
  T, result_t,
  mpl::enable_if_t<
    std::is_same<
      result_t,
      decltype(
         std::declval<T const>().createResidual()
         )
      >::value
    >
  > : std::true_type{};

}}} 
#endif
