
#ifndef ROM_HAS_CONST_CREATE_apply_TD_JAC_OBJECT_result_METHOD_RETURN_RESULT_HPP_
#define ROM_HAS_CONST_CREATE_apply_TD_JAC_OBJECT_result_METHOD_RETURN_RESULT_HPP_

namespace pressio{ namespace rom{ namespace predicates {

template <
  typename T, typename operand_t, typename result_t,
  typename = void
  >
struct has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result
  : std::false_type{};


template <typename T, typename operand_t, typename result_t>
struct has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result<
  T, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    mpl::is_same<
      result_t,
      decltype(
	       std::declval<T const>().createApplyDiscreteTimeJacobianResult
          (
            std::declval<operand_t const & >()
          )
	       )
      >::value
    >
  > : std::true_type{};


}}} 
#endif
