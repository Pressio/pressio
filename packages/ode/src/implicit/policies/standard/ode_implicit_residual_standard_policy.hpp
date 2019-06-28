
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace rompp{ namespace ode{ namespace policy{

/*
 * state and residual_types are algebra wrappers
 * both are wrappers from algebra
 */
template<
  typename state_type,
  typename model_type,
  typename residual_type
  >
class ImplicitResidualStandardPolicy<
  state_type, model_type, residual_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value and
    algebra::meta::is_wrapper<state_type>::value and
    algebra::meta::is_wrapper<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicy<state_type, model_type, residual_type>>
{

  using this_t = ImplicitResidualStandardPolicy<state_type, model_type, residual_type>;
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
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const{

    model.residual(*y.data(), *R.data(), t);
    ::rompp::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
  }

  template <
    ode::ImplicitEnum method, int n, typename scalar_type
    >
  residual_type operator()(const state_type & y,
  			   const std::array<state_type, n> & oldYs,
  			   const model_type & model,
  			   scalar_type t,
  			   scalar_type dt)const {

    auto nR = model.residual(*y.data(), t);
    residual_type R(nR);
    ::rompp::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
    return R;
  }

};//end class

}}}//end namespace rompp::ode::policy
#endif
