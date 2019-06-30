
#ifndef CONTAINERS_META_HAS_STATIC_METHOD_SCALE_HPP_
#define CONTAINERS_META_HAS_STATIC_METHOD_SCALE_HPP_

#include <utility>

namespace rompp{ namespace containers{ namespace meta {

template <
  typename T,
  typename arg_t,
  typename scalar_t,
  typename = void
  >
struct has_static_method_scale : std::false_type{};

template <
  typename T,
  typename arg_t,
  typename sc_t
  >
struct has_static_method_scale<
  T, arg_t, sc_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       T::scale
       (
	std::declval< arg_t & >(),
	std::declval< const sc_t >()
	)
       )
      >::value
    >
  > : std::true_type{};

}}} //rompp::containers::meta
#endif
