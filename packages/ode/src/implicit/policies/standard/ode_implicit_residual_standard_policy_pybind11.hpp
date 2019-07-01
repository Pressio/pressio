
#ifdef HAVE_PYBIND11
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_PYBIND11_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_PYBIND11_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace rompp{ namespace ode{ namespace policy{

template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ImplicitResidualStandardPolicyPybind11<
  state_type, system_type, residual_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value and
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_cstyle_array_pybind11<state_type>::value and
    containers::meta::is_cstyle_array_pybind11<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicyPybind11<state_type, system_type, residual_type>>
{

  using this_t = ImplicitResidualStandardPolicyPybind11<state_type, system_type, residual_type>;
  friend ImplicitResidualPolicyBase<this_t>;

public:
  ImplicitResidualStandardPolicyPybind11() = default;
  ~ImplicitResidualStandardPolicyPybind11() = default;

public:
  template <
    ode::ImplicitEnum method, int n, typename scalar_type
  >
  void operator()(const state_type & y,
		  residual_type & R,
		  const std::array<state_type, n> & oldYs,
		  const system_type & model,
		  scalar_type t,
		  scalar_type dt) const{
    
    throw std::runtime_error("ImplicitResidualStandardPolicyPybind11 missing");

    // printf("C++ R address: %p\n", R.data());
    // model.attr("residual2")(y, R, t);
    // ::rompp::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
  }

  template <
    ode::ImplicitEnum method, int n, typename scalar_type
    >
  residual_type operator()(const state_type & y,
  			   const std::array<state_type, n> & oldYs,
  			   const system_type & model,
  			   scalar_type t,
  			   scalar_type dt) const {

    throw std::runtime_error("ImplicitResidualStandardPolicyPybind11 missing");

    residual_type nR;// = model.attr("residual1")(y, t);
    // ::rompp::ode::impl::time_discrete_residual<method, n>(y, nR, oldYs, dt);
    return nR;
  }
};//end class

}}}//end namespace rompp::ode::policy
#endif
#endif
