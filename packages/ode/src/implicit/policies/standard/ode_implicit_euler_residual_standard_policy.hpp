
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../../ode_forward_declarations.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"
#include "../../../ode_basic_meta.hpp"
#include "../../ode_implicit_constants.hpp"

namespace rompp{ namespace ode{ namespace policy{

template<typename state_type,
	 typename model_type,
	 typename residual_type>
class ImplicitEulerResidualStandardPolicy<
  state_type, model_type, residual_type,
  core::meta::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
            ImplicitEulerResidualStandardPolicy<state_type, model_type, residual_type>,
            ::rompp::ode::coeffs::bdf1_numAuxStates_,
            ::rompp::ode::coeffs::bdf1_numAuxRHS_>
{

  using this_t = ImplicitEulerResidualStandardPolicy<state_type, model_type, residual_type>;
  using arr_states_t = std::array<state_type, ::rompp::ode::coeffs::bdf1_numAuxStates_>;

public:
  ImplicitEulerResidualStandardPolicy() = default;
  ~ImplicitEulerResidualStandardPolicy() = default;

private:
  using scalar_type =
    typename core::details::traits<state_type>::scalar_t;

public:
  //----------------------------------------------------------------
  void operator()(const state_type & y,
		  residual_type & R,
		  const arr_states_t & oldYs,
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const{
    R.setZero();
    model.residual(*y.data(), *R.data(), t);
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }
  //----------------------------------------------------------------

  residual_type operator()(const state_type & y,
			   const arr_states_t & oldYs,
			   const model_type & model,
			   scalar_type t,
			   scalar_type dt)const {

    auto nR = model.residual(*y.data(), t);
    residual_type R(nR);
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
    return R;
  }
  //----------------------------------------------------------------

private:
  friend ImplicitResidualPolicyBase<
  ImplicitEulerResidualStandardPolicy<state_type, model_type, residual_type>,
  ::rompp::ode::coeffs::bdf1_numAuxStates_, ::rompp::ode::coeffs::bdf1_numAuxRHS_>;
};//end class

}}}//end namespace rompp::ode::policy
#endif
