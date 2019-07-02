
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"

namespace pressio{ namespace ode{ namespace policy{

/*
 * state and jacobian types are containers wrappers
 * both are wrappers from containers
 */
template<typename state_type,
	 typename system_type,
	 typename jacobian_type>
class ImplicitJacobianStandardPolicy<
  state_type, system_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
    state_type, system_type, jacobian_type> >{

  using this_t = ImplicitJacobianStandardPolicy<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicy() = default;
  ~ImplicitJacobianStandardPolicy() = default;

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
    model.jacobian( *y.data(), *J.data(), t);
    ::pressio::ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <
    ode::ImplicitEnum method, typename scalar_t
    >
  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   scalar_t t,
  			   scalar_t dt)const
  {
    auto nJJ = model.jacobian(*y.data(), t);
    jacobian_type JJ(nJJ);
    ::pressio::ode::impl::time_discrete_jacobian<method>(JJ, dt);
    return JJ;
  }

};//end class



#ifdef HAVE_PYBIND11
/*
 * state_type and jacobian_type = pybind11::array_t
 */
template<typename state_type,
	 typename system_type,
	 typename jacobian_type>
class ImplicitJacobianStandardPolicy<
  state_type, system_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_array_pybind11<state_type>::value and
    containers::meta::is_array_pybind11<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
    state_type, system_type, jacobian_type> >{

  using this_t = ImplicitJacobianStandardPolicy<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicy() = default;
  ~ImplicitJacobianStandardPolicy() = default;

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
    throw std::runtime_error(" ImplicitJacobianStandardPolicy for pybind11 missing ");
    // model.jacobian( *y.data(), *J.data(), t);
    // ::pressio::ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <
    ode::ImplicitEnum method, typename scalar_t
    >
  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   scalar_t t,
  			   scalar_t dt)const
  {
    auto nJJ = model.attr("jacobian1")(y, t);
    throw std::runtime_error(" ImplicitJacobianStandardPolicy for pybind11 missing ");

    // auto nJJ = model.jacobian(*y.data(), t);
    // jacobian_type JJ(nJJ);
    // ::pressio::ode::impl::time_discrete_jacobian<method>(JJ, dt);
    return nJJ;
  }
};//end class
#endif


}}}//end namespace pressio::ode::policy
#endif
