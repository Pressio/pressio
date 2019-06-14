
#ifndef ODE_IMPLICIT_AUX_DATA_HPP_
#define ODE_IMPLICIT_AUX_DATA_HPP_

#include "ode_ConfigDefs.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace rompp{ namespace ode{ namespace impl{

template<
  typename model_type,
  typename scalar_type
  >
struct ImplicitOdeAuxData<
  model_type, scalar_type
#ifdef HAVE_PYBIND11
  , mpl::enable_if_t<
      ::rompp::mpl::not_same<model_type, pybind11::object >::value
      >
#endif
  >
{
  const model_type & model_;
  scalar_type t_  = {};
  scalar_type dt_ = {};

  ImplicitOdeAuxData(const model_type & model)
    : model_(model){}

  ~ImplicitOdeAuxData() = default;
};


#ifdef HAVE_PYBIND11
template<
  typename model_type,
  typename scalar_type
  >
struct ImplicitOdeAuxData<
  model_type, scalar_type,
  mpl::enable_if_t<
    ::rompp::mpl::is_same<model_type, pybind11::object >::value
    >
  >
{
  model_type model_;
  scalar_type t_  = {};
  scalar_type dt_ = {};

  ImplicitOdeAuxData(const model_type & model)
    : model_(model){}

  ~ImplicitOdeAuxData() = default;
};
#endif

}}}//end namespace rompp::ode::impl
#endif
