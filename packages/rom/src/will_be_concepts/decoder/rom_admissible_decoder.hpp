
#ifndef ROM_ADMISSIBLE_DECODER_HPP_
#define ROM_ADMISSIBLE_DECODER_HPP_

namespace pressio{ namespace rom{ namespace concepts {

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
    ::pressio::rom::predicates::has_const_get_reference_to_jacobian<T, typename T::jacobian_type>::value and
    ::pressio::rom::predicates::has_const_apply_mapping_accept_operand_result_return_void<T, operand_t, result_t>::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::concepts
#endif
