
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../../ode_forward_declarations.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"
#include "../../../ode_basic_meta.hpp"

namespace rompp{ namespace ode{ namespace policy{

template<typename state_type,
	 typename model_type,
	 typename residual_type>
class ImplicitResidualStandardPolicy<
  state_type, model_type, residual_type,
  core::meta::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicy<state_type, model_type, residual_type>>
{

  using this_t = ImplicitResidualStandardPolicy<state_type, model_type, residual_type>;
  using scalar_type = typename core::details::traits<state_type>::scalar_t;

public:
  ImplicitResidualStandardPolicy() = default;
  ~ImplicitResidualStandardPolicy() = default;

public:

  template <::rompp::ode::ImplicitEnum method,
	     int numAuxStates>
  void operator()(const state_type & y,
		  residual_type & R,
		  const std::array<state_type, numAuxStates> & oldYs,
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const{

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
    ode::impl::implicit_time_discrete_residual<method, numAuxStates>(y, oldYs, R, dt);
  }
  //----------------------------------------------------------------

  template <::rompp::ode::ImplicitEnum method,
	     int numAuxStates>
  residual_type operator()(const state_type & y,
			   const std::array<state_type, numAuxStates> & oldYs,
			   const model_type & model,
			   scalar_type t,
			   scalar_type dt)const {

    auto nR = model.residual(*y.data(), t);
    residual_type R(nR);
    ode::impl::implicit_time_discrete_residual<method, numAuxStates>(y, oldYs, R, dt);
    return R;
  }
  //----------------------------------------------------------------

private:
  friend ImplicitResidualPolicyBase<this_t>;

};//end class

}}}//end namespace rompp::ode::policy
#endif
