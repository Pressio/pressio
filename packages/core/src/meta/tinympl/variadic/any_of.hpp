
#ifndef TINYMPL_VARIADIC_ANY_OF_HPP
#define TINYMPL_VARIADIC_ANY_OF_HPP

#include <type_traits>

namespace tinympl { namespace variadic {

/**
 * \ingroup VarNonModAlgs
 * \class any_of
 * \brief Determines whether any of the elements in the sequence satisfy the given predicate
 * \param Predicate the predicate, `Predicate<T>::type::value` must be convertible to `bool`
 * \param Args... the input sequence
 * \return `any_of<...>::type` is a `std::integral_constant<bool,v>` where `v`
is true iff at least one element in the sequence satisfy the predicate `Predicate`
 * \sa tinympl::any_of
 */
template< template<class... T> class Predicate, class ... Args>
struct any_of;

template< template<class... T> class Predicate, class Head, class... Args>
struct any_of<Predicate, Head, Args...>
  : std::conditional<
  Predicate<Head>::type::value,
  std::integral_constant<bool, true>,
  typename any_of<Predicate, Args...>::type
  >::type
{};

template< template<class ... T> class Predicate>
struct any_of<Predicate>
  : std::integral_constant<bool, false>{};


}} // namespace tinympl::variadic

#endif // TINYMPL_VARIADIC_ANY_OF_HPP
