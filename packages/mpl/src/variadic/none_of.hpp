
#ifndef ROMPP_MPL_VARIADIC_NONE_OF_HPP
#define ROMPP_MPL_VARIADIC_NONE_OF_HPP

#include <type_traits>

namespace rompp{ namespace mpl{ namespace variadic {

/**
 * \class none_of
 * \brief Determines whether none of the elements in the sequence satisfy the given predicate
 * \param Predicate the predicate, `Predicate<T>::type::value` must be convertible to bool
 * \param Args... the input sequence
 * \return `none_of<...>::type` is a `std::integral_constant<bool,v>` where `v`
is true iff none of the elements in the sequence satisfy the predicate `Predicate`
 * \sa tinympl::none_of
 */
template< template<class ... T> class Predicate, class ... Args>
struct none_of;

template< template<class ... T> class Predicate, class Head, class ... Args>
struct none_of<Predicate, Head, Args...>
  : std::conditional <
  Predicate<Head>::type::value,
  std::integral_constant<bool, false>,
  typename none_of<Predicate, Args...>::type
  >::type
{};

template< template<class ... T> class Predicate>
struct none_of<Predicate>
  : std::integral_constant<bool, true>{};


}}} // namespace rompp::mpl::variadic

#endif 
