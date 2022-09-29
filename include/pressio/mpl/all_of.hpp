
#ifndef MPL_ALL_OF_HPP_
#define MPL_ALL_OF_HPP_

#include "./variadic/all_of.hpp"

namespace pressio{ namespace mpl{

template< template<class ... T> class F, class ... Args>
struct all_of : variadic::all_of<F, Args...>{};

template< template<class... T> class F, class ... Args>
using all_of_t = typename all_of<F, Args...>::type;

}}
#endif  // MPL_ALL_OF_HPP_
