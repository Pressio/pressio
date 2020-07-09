
#ifndef rom_has_const_apply_mapping_accept_operand_result_return_void_HPP_
#define rom_has_const_apply_mapping_accept_operand_result_return_void_HPP_

namespace pressio{ namespace rom{ namespace predicates {

template <typename T, typename operand_t, typename result_t, typename = void>
struct has_const_apply_mapping_accept_operand_result_return_void : std::false_type{};

template <typename T, typename operand_t, typename result_t>
struct has_const_apply_mapping_accept_operand_result_return_void<
  T, operand_t, result_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const &>().applyMapping
       (
	     std::declval<operand_t const &>(), 
         std::declval<result_t &>()
	    )
       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::predicates
#endif
