

#ifndef ode_has_const_create_velocity_method_return_result_HPP_
#define ode_has_const_create_velocity_method_return_result_HPP_

namespace pressio{ namespace ode{ namespace meta {

template <typename T, typename velo_type, typename = void>
struct has_const_create_velocity_method_return_result
  : std::false_type{};

template <typename T, typename velo_type>
struct has_const_create_velocity_method_return_result<
  T, velo_type,
  ::pressio::mpl::enable_if_t<
    !std::is_void<velo_type>::value and
    mpl::is_same<
      velo_type,
      decltype(
	       std::declval<T const>().createVelocity()
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
