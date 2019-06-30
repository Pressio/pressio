
#ifndef ODE_EXPLICIT_AUX_DATA_HPP_
#define ODE_EXPLICIT_AUX_DATA_HPP_

#include "ode_ConfigDefs.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace rompp{ namespace ode{ namespace impl{

template<
  typename model_type,
  typename residual_policy_type
  >
class ExplicitOdeAuxData<
  model_type, residual_policy_type
#ifdef HAVE_PYBIND11
  , mpl::enable_if_t<
      ::rompp::mpl::not_same<model_type, pybind11::object >::value
      >
#endif
  >{
protected:
  ExplicitOdeAuxData(const model_type & mod,
		     const residual_policy_type & rpolo)
    : model_(mod), residual_obj_(&rpolo){}

  ~ExplicitOdeAuxData() = default;

  const model_type & model_;
  const residual_policy_type * residual_obj_ = nullptr;
};



/* for some reason to be determined, when we deal with
 * python objects, we need to pass by copy
 */
#ifdef HAVE_PYBIND11
template<
  typename model_type,
  typename residual_policy_type
  >
class ExplicitOdeAuxData<
  model_type, residual_policy_type,
  mpl::enable_if_t<
    ::rompp::mpl::is_same<model_type, pybind11::object >::value
    >
  >{
protected:
  ExplicitOdeAuxData(const model_type mod,
		     const residual_policy_type & rpolo)
    : model_(mod), residual_obj_(&rpolo){}

  ~ExplicitOdeAuxData() = default;

  model_type model_;
  const residual_policy_type * residual_obj_ = nullptr;
};
#endif

}}}//end namespace rompp::ode::impl
#endif
