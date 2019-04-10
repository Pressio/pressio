
#ifndef TINYMPL_VARIADIC_SIZE_HPP
#define TINYMPL_VARIADIC_SIZE_HPP

#include <type_traits>
#include <cstddef>

namespace tinympl { namespace variadic {

/**
 * \ingroup VariadicBasic
 * \class size
 * \brief Compute the size of a variadic template
 * \return `size<Args...>::value` is equivalent to `sizeof ... (Args)`
 */
template<class ... Args> struct size;

template<class ... Args> using size_t = typename size<Args...>::type;

template<class ... Args> struct size : std::integral_constant<std::size_t,sizeof...(Args)> {};

}} // namespace tinympl::variadic

#endif // TINYMPL_VARIADIC_SIZE_HPP
