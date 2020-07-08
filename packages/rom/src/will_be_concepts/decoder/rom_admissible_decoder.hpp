
#ifndef ROM_ADMISSIBLE_DECODER_HPP_
#define ROM_ADMISSIBLE_DECODER_HPP_

namespace pressio{ namespace rom{ namespace meta {

template <typename T, typename jacobian_t, typename = void>
struct has_get_reference_to_jacobian : std::false_type{};

template <typename T, typename jacobian_t>
struct has_get_reference_to_jacobian<
  T, jacobian_t,
  mpl::enable_if_t<
    mpl::is_same<
      decltype( std::declval<T const &>().getReferenceToJacobian() ),
      const jacobian_t &
      >::value
    >
  > : std::true_type{};
//--------------------------------------------------------------


template <typename T, typename arg1_t, typename arg2_t, typename = void>
struct has_apply_mapping_two_args : std::false_type{};

template <typename T, typename arg1_t, typename arg2_t>
struct has_apply_mapping_two_args<
  T, arg1_t, arg2_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const &>().applyMapping
       (
	     std::declval<arg1_t const &>(), 
       std::declval<arg2_t &>()
	     )
       )
      >::value
    >
  > : std::true_type{};
//--------------------------------------------------------------


/*
 * A type is a legitimate decoder for LSPG if:
 *
 * - has a jacobian_typedef
 * - has a getReferenceToJacobian
 * - template applyMapping(operand_t, result_t)
 *
*/
template<
  typename T,
  typename operand_t,
  typename result_t,
  typename enable = void
  >
struct admissible_decoder : std::false_type{};

template<
  typename T,
  typename operand_t,
  typename result_t
  >
struct admissible_decoder<
  T, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::predicates::has_jacobian_typedef<T>::value and
    has_get_reference_to_jacobian<T, typename T::jacobian_type>::value and
    has_apply_mapping_two_args<T, operand_t, result_t>::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::meta
#endif
