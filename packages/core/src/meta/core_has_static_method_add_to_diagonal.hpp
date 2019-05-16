
#ifndef CORE_META_HAS_STATIC_METHOD_ADD_TO_DIAGONAL_HPP_
#define CORE_META_HAS_STATIC_METHOD_ADD_TO_DIAGONAL_HPP_

#include <utility>

namespace rompp{ namespace core{ namespace meta {

template <
  typename T,
  typename arg_t,
  typename scalar_t,
  typename = void
  >
struct has_static_method_add_to_diagonal : std::false_type{};

template <
  typename T,
  typename arg_t,
  typename sc_t
  >
struct has_static_method_add_to_diagonal<
  T, arg_t, sc_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       T::add_to_diagonal
       (
	std::declval< arg_t & >(),
	std::declval< const sc_t >()
	)
       )
      >::value
    >
  > : std::true_type{};

}}} //rompp::core::meta
#endif
