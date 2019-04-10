#ifndef ROMPP_MPL_VARIADIC_AT_HPP
#define ROMPP_MPL_VARIADIC_AT_HPP

#include <cstddef>
#include "../identity.hpp"
#include <tuple>

namespace rompp{ namespace mpl{ namespace variadic {

/**
 * \ingroup VarBasics
 * \class at
 * \brief Extract the i-th element of a variadic template
 * \param i The index to extract
 */
template <std::size_t i, typename... Args> struct at;

template <std::size_t i, typename... Args>
struct at
  : identity<
  typename std::tuple_element<i, std::tuple<Args...>>::type
  >{};

template<std::size_t i, class... Args>
using at_t = typename at<i, Args...>::type;
//----------------------------------------------------------------

template<typename Default, std::size_t i, typename ... Args>
struct at_or
  : std::conditional<
  i < sizeof ... (Args),
      at<i, Args...>,
      identity<Default>
      >::type
  {};

template<typename Default, std::size_t i, typename ... Args>
using at_or_t = typename at_or<Default, i, Args...>::type;

}}} // namespace 

#endif // ROMPP_MPL_VARIADIC_AT_HPP




// template<std::size_t i, typename ... Args> struct at;

// template<std::size_t i, typename Head, typename ... Tail>
// struct at<i, Head, Tail...>{
//   static_assert(i < sizeof ... (Tail) + 1, "Index out of range");
//   using type = typename at<i-1,Tail...>::type;
// };

// template<typename Head, typename ... Tail>
// struct at<0,Head,Tail...>{
//   using type = Head;
// };

// template<std::size_t i, typename ... Args>
// using at_t = typename at<i, Args...>::type;
