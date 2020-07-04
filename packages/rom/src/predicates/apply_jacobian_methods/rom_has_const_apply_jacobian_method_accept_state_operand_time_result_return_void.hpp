
#ifndef rom_has_const_apply_jacobian_method_accept_state_operand_time_result_return_void_HPP_
#define rom_has_const_apply_jacobian_method_accept_state_operand_time_result_return_void_HPP_

namespace pressio{ namespace rom{ namespace meta {

template <
  typename T,
  typename state_t,
  typename operand_t,
  typename time_type,
  typename result_t,
  typename = void
  >
struct has_const_apply_jacobian_method_accept_state_operand_time_result_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename operand_t,
  typename time_type,
  typename result_t
  >
struct has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
  T, state_t, operand_t, time_type, result_t,
  ::pressio::mpl::void_t<
  decltype(
	   std::declval<T const>().applyJacobian(
					    std::declval<state_t const&>(),
              std::declval<operand_t const&>(),
					    std::declval<time_type const &>(),
					    std::declval<result_t &>()
					    )
	   )
    >
  >: std::true_type{};

}}} 
#endif
