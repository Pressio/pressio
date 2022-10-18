
#ifndef MPL_AT_HPP_
#define MPL_AT_HPP_

#include "./variadic/at.hpp"

namespace pressio{ namespace mpl{

template <std::size_t i, typename... Args>
struct at : variadic::at<i, Args...>{};

template<std::size_t i, class... Args>
using at_t = typename at<i, Args...>::type;

template<typename Default, std::size_t i, typename ... Args>
struct at_or : variadic::at_or<Default, i, Args...>{};

template<typename Default, std::size_t i, typename ... Args>
using at_or_t = typename at_or<Default, i, Args...>::type;

}}

#endif  // MPL_AT_HPP_
