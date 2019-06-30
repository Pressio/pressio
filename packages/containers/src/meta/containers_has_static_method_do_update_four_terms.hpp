
#ifndef CONTAINERS_HAS_STATIC_METHOD_DO_UPDATE_FOUR_TERMS_HPP_
#define CONTAINERS_HAS_STATIC_METHOD_DO_UPDATE_FOUR_TERMS_HPP_

#include <utility>

namespace rompp{ namespace containers{ namespace meta {

/*
  detect if type T has a static method of the form:

  static void do_update(T1 &      , scalar_type,
			const T2 &, scalar_type,
			const T3 &, scalar_type,
			const T4 &, scalar_type,
			const T5 &, scalar_type)
 */

template <
  typename T,
  typename scalar_t,
  typename T1,
  typename T2,
  typename T3,
  typename T4,
  typename T5,
  typename = void
  >
struct has_static_method_do_update_four_terms : std::false_type{};

template <
  typename T,
  typename sc_t,
  typename T1,
  typename T2,
  typename T3,
  typename T4,
  typename T5
  >
struct has_static_method_do_update_four_terms<
  T, sc_t, T1, T2, T3, T4, T5,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       T::do_update
       (
	std::declval< T1 & >(),
	std::declval<const sc_t>(),
	std::declval<const T2 &>(),
	std::declval<const sc_t>(),
	std::declval<const T3 &>(),
	std::declval<const sc_t>(),
	std::declval<const T4 &>(),
	std::declval<const sc_t>(),
	std::declval<const T5 &>(),
	std::declval<const sc_t>()
	)
       )
      >::value
    >
  > : std::true_type{};

}}} // namespace rompp::containers::meta
#endif
