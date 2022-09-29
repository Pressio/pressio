
#ifndef MPL_NONE_OF_HPP_
#define MPL_NONE_OF_HPP_

#include "./variadic/none_of.hpp"

namespace pressio{ namespace mpl{

template< template<class ... T> class F, class ... Args>
struct none_of : variadic::none_of<F, Args...>{};

template< template<class... T> class F, class ... Args>
using none_of_t = typename none_of<F, Args...>::type;

}}
#endif  // MPL_NONE_OF_HPP_
