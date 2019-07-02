
#ifndef PRESSIO_MPL_VARIADIC_SIZE_HPP
#define PRESSIO_MPL_VARIADIC_SIZE_HPP

#include <type_traits>
#include <cstddef>

namespace pressio{ namespace mpl{ namespace variadic {

/**
 * \class size
 * \brief Compute the size of a variadic template
 * \return `size<Args...>::value` is equivalent to `sizeof ... (Args)`
 */
template<class ... Args> 
struct size;

template<class ... Args> 
using size_t = typename size<Args...>::type;

template<class ... Args> 
struct size : std::integral_constant<std::size_t,sizeof...(Args)>
{};

}}} // namespace pressio::mpl::variadic

#endif 
