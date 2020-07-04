
#ifndef PRESSIO_SOLVERS_NEWTON_RAPHSON_HPP_
#define PRESSIO_SOLVERS_NEWTON_RAPHSON_HPP_

#include "./impl/solvers_nonlinear_compose.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{

template<
  typename system_t,
  template<typename...> class update,
  template<typename...> class looper,
  typename ... Args
  >
using composeNewtonRaphson = impl::compose<system_t, impl::NewtonRaphson, update, looper, Args...>;

template<
  typename system_t,
  template<typename...> class update,
  template<typename...> class looper,
  typename ... Args
  >
using composeNewtonRaphson_t = typename composeNewtonRaphson<system_t, update, looper, Args...>::type;

}}}
#endif
