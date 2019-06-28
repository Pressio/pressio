
#ifndef ALGEBRA_HAS_STATIC_METHOD_DO_UPDATE_ONE_TERM_HPP_
#define ALGEBRA_HAS_STATIC_METHOD_DO_UPDATE_ONE_TERM_HPP_

#include <utility>

namespace rompp{ namespace algebra{ namespace meta {

/*
 * detect if type T has a static method of the form:
 *
 * static void do_update(T1 & , scalar_type, const T2 &, scalar_type)
 */

template <typename T,
	  typename scalar_t,
	  typename T1,
	  typename T2,
	  typename = void>
struct has_static_method_do_update_one_term : std::false_type{};

template <typename T,
	  typename sc_t,
	  typename T1,
	  typename T2 >
struct has_static_method_do_update_one_term<
  T, sc_t, T1, T2,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       T::do_update
       (
	std::declval< T1 & >(),
	std::declval<const sc_t>(),
	std::declval<const T2 &>(),
	std::declval<const sc_t>()
	)
       )
      >::value
    and
    std::is_void<
      decltype
      (
       T::do_update
       (
	std::declval< T1 & >(),
	std::declval<const T2 &>(),
	std::declval<const sc_t>()
	)
       )
      >::value
    >
  > : std::true_type{};

// template <typename T,
// 	  typename sc_t,
// 	  typename T1,
// 	  typename T2 >
// struct has_static_method_do_update_one_term<
//   T, sc_t, T1, T2,
//   mpl::enable_if_t<
//     std::is_void<
//       decltype
//       (
//        T::do_update
//        (
// 	std::declval< T1 & >(),
// 	std::declval<const T2 &>(),
// 	std::declval<const sc_t>()
// 	)
//        )
//       >::value
//     >
//   > : std::true_type{};

}}} // namespace rompp::algebra::meta
#endif
