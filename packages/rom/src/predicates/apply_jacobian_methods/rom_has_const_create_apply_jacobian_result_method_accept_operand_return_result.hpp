
#ifndef rom_has_const_create_apply_jacobian_result_method_accept_operand_return_result_HPP_
#define rom_has_const_create_apply_jacobian_result_method_accept_operand_return_result_HPP_

namespace pressio{ namespace rom{ namespace predicates {
  
template <
  typename T,
  typename operand_t,
  typename result_t,
  typename = void
  >
struct has_const_create_apply_jacobian_result_method_accept_operand_return_result
  : std::false_type{};

template <
  typename T,
  typename operand_t,
  typename result_t
  >
struct has_const_create_apply_jacobian_result_method_accept_operand_return_result<
  T, operand_t, result_t,
  mpl::enable_if_t<
    std::is_same<
      result_t,
      decltype(
         std::declval<T const>().createApplyJacobianResult
          (
          std::declval<operand_t const &>()
          )
         )
      >::value
    >
  > : std::true_type{};

}}} 
#endif
