
#ifdef HAVE_PYBIND11
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_PYBIND11_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_PYBIND11_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"

namespace rompp{ namespace ode{ namespace policy{

template<typename state_type,
	 typename system_type,
	 typename jacobian_type>
class ImplicitJacobianStandardPolicyPybind11<
  state_type, system_type, jacobian_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_cstyle_array_pybind11<state_type>::value and
    containers::meta::is_cstyle_array_pybind11<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicyPybind11<
    state_type, system_type, jacobian_type> >{

  using this_t = ImplicitJacobianStandardPolicyPybind11<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicyPybind11() = default;
  ~ImplicitJacobianStandardPolicyPybind11() = default;

public:
  template <
    ode::ImplicitEnum method, typename scalar_t
  >
  void operator()(const state_type & y,
		  jacobian_type & J,
		  const system_type & model,
		  scalar_t t,
		  scalar_t dt)const
  {
    throw std::runtime_error("ImplicitJacobianStandardPolicyPybind11 is missing");    
    // model.attr("jacobian2")(y, J, t);
    // ::rompp::ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <
    ode::ImplicitEnum method, typename scalar_t
    >
  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   scalar_t t,
  			   scalar_t dt)const
  {
    throw std::runtime_error("ImplicitJacobianStandardPolicyPybind11 is missing");    
    jacobian_type nJJ;// = model.attr("jacobian1")(y, t);
    // ::rompp::ode::impl::time_discrete_jacobian<method>(nJJ, dt);
    return nJJ;
  }
};//end class

}}}//end namespace rompp::ode::policy
#endif
#endif
