
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../ode_forward_declarations.hpp"
#include "ode_explicit_residual_policy_base.hpp"

namespace rompp{ namespace ode{ namespace policy{

  
// state and residual have same type and are
// wrappers from core
template<typename state_type,
	 typename model_type>
class ExplicitResidualStandardPolicy<
  state_type,model_type, state_type,
    core::meta::enable_if_t<
      // enable when state and residual are vector wrappers
      // what about the case when they are multivector wrappers?
      // think if this works right away
      core::meta::is_core_vector_wrapper<state_type>::value
      >
  >
  : public ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<state_type, model_type> >{

public:
  ExplicitResidualStandardPolicy() = default;
  ~ExplicitResidualStandardPolicy() = default;  
  
  using scalar_type =
    typename core::details::traits<state_type>::scalar_t;
  
  void operator()(const state_type & y,
		  state_type & R,
		  model_type & model,
		  scalar_type t){
    if (R.empty())
      R.matchLayoutWith(y);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
  }
  //----------------------------------------------
  
private:
  friend ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<
    state_type, model_type> >;

};//end class

}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 
