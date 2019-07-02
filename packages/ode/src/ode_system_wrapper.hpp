
#ifndef ODE_SYSTEM_WRAPPER_HPP_
#define ODE_SYSTEM_WRAPPER_HPP_

#include "ode_ConfigDefs.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace pressio{ namespace ode{ namespace impl{

template<typename model_type>
struct OdeSystemWrapper<
  model_type
#ifdef HAVE_PYBIND11
  , mpl::enable_if_t<
      ::pressio::mpl::not_same<model_type, pybind11::object >::value
      >
#endif
  >
{
  OdeSystemWrapper(const model_type & system)
    : data_(system){}

  OdeSystemWrapper() = delete;
  ~OdeSystemWrapper() = default;

  const model_type & get() const{
    return data_;
  }

private:
  const model_type & data_;
};


/* for some reason to be determined, when we deal with
 * python objects, we need to pass by copy
 */
#ifdef HAVE_PYBIND11
template<typename model_type>
class OdeSystemWrapper<
  model_type,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<model_type, pybind11::object >::value
    >
  >
{
  OdeSystemWrapper(const model_type & system)
    : data_(system){}

  OdeSystemWrapper() = delete;
  ~OdeSystemWrapper() = default;

  const model_type & get() const{
    return data_;
  }

private:
    model_type data_;
};
#endif

}}}//end namespace pressio::ode::impl
#endif
