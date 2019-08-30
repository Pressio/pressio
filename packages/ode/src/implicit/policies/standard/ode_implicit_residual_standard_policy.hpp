
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace pressio{ namespace ode{ namespace policy{

/*
 * state and residual_types are containers wrappers
 * both are wrappers from containers
 */
template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ImplicitResidualStandardPolicy<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicy<state_type, system_type, residual_type>>
{

  using this_t = ImplicitResidualStandardPolicy<state_type, system_type, residual_type>;
  friend ImplicitResidualPolicyBase<this_t>;

public:
  ImplicitResidualStandardPolicy() = default;
  ~ImplicitResidualStandardPolicy() = default;

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

    model.velocity(*y.data(), t, *R.data());
    ::pressio::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
  }

  template <
    ode::ImplicitEnum method, int n, typename scalar_type
    >
  residual_type operator()(const state_type & y,
  			   const std::array<state_type, n> & oldYs,
  			   const system_type & model,
  			   scalar_type t,
  			   scalar_type dt)const {

    auto nR = model.velocity(*y.data(), t);
    residual_type R(nR);
    ::pressio::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
    return R;
  }

};//end class

}}}//end namespace pressio::ode::policy
#endif
