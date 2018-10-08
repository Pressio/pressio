
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_BDF2_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_BDF2_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../../ode_forward_declarations.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace rompp{ namespace ode{ namespace policy{
  
//----------------------------------------------------------------
//   default when state_type = residual_type
//----------------------------------------------------------------
template<typename state_type,
	 typename model_type>
class ImplicitBDF2ResidualStandardPolicy<
  state_type, model_type, state_type,
  core::meta::enable_if_t<
    core::meta::is_core_vector_wrapper<state_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitBDF2ResidualStandardPolicy<state_type,model_type>,2,0>{
  
public:
  ImplicitBDF2ResidualStandardPolicy() = default;
  ~ImplicitBDF2ResidualStandardPolicy() = default;  

private:
  using scalar_type =
    typename core::details::traits<state_type>::scalar_t;

public:  
  //----------------------------------------------------------------
  void operator()(const state_type & y,
		  state_type & R,
		  const std::array<state_type, 2> & oldYs,
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const{    
    if (R.empty())
      R.matchLayoutWith(y);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
    // do time discrete residual
    ode::impl::implicit_bdf2_time_discrete_residual(y,
						    oldYs[1],
						    oldYs[0],
						    R, dt);
  }
  //----------------------------------------------------------------

  state_type operator()(const state_type & y, 
			const std::array<state_type, 2> & oldYs,
			const model_type & model,
			scalar_type t,
			scalar_type dt)const {
    
    auto nR = model.residual(*y.data(), t);
    state_type R(nR);
    ode::impl::implicit_bdf2_time_discrete_residual(y,
						    oldYs[1],
						    oldYs[0],
						    R, dt);
    return R;
  }
  //----------------------------------------------------------------
  
private:
  friend ImplicitResidualPolicyBase<
    ImplicitBDF2ResidualStandardPolicy<state_type, model_type,
					state_type>, 2,0>;
};//end class

      
}}}//end namespace rompp::ode::policy
#endif 
