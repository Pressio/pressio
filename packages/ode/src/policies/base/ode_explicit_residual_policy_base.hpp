
#ifndef ODE_POLICIES_BASE_EXPLICIT_RESIDUAL_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_EXPLICIT_RESIDUAL_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
    
template <typename derived_t>
class ExplicitResidualPolicyBase
  : private core::details::CrtpBase<ExplicitResidualPolicyBase<derived_t>>
{
public:
  template <typename state_type,
	    typename residual_type,
	    typename model_type,
	    typename time_type>
  void compute(const state_type & y,
	       residual_type & R,
	       model_type & model,
	       time_type t){
    this->underlying().computeImpl(y, R, model, t);
  }

private:
  friend derived_t;
  friend core::details::CrtpBase<ExplicitResidualPolicyBase<derived_t>>;

  ExplicitResidualPolicyBase() = default;
  ~ExplicitResidualPolicyBase() = default;
  

};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
