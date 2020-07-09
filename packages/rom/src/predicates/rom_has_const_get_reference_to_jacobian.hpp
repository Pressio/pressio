
#ifndef rom_has_const_get_reference_to_jacobian_HPP_
#define rom_has_const_get_reference_to_jacobian_HPP_

namespace pressio{ namespace rom{ namespace predicates {

template <typename T, typename jacobian_t, typename = void>
struct has_const_get_reference_to_jacobian : std::false_type{};

template <typename T, typename jacobian_t>
struct has_const_get_reference_to_jacobian<
  T, jacobian_t,
  mpl::enable_if_t<
    mpl::is_same<
      decltype( std::declval<T const &>().getReferenceToJacobian() ),
      const jacobian_t &
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::predicates
#endif
