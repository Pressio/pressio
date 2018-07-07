
#ifndef ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../../base/ode_residual_policy_base.hpp"
#include "../../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, 
	typename residual_type,
	typename model_type, 
	typename time_type, typename sizer_type>
class implicitEulerStandardResidual
  : public residualPolicyBase<implicitEulerStandardResidual,
					   state_type, residual_type,
					   model_type, time_type, sizer_type>,
  public advanceFullStatePolicyBase<implicitEulerStandardResidual,
				    state_type, residual_type,
				    model_type, time_type, sizer_type>
{
public:
  implicitEulerStandardResidual() = default;
  ~implicitEulerStandardResidual() = default;  

private:
  // enable if using types from core package
  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreVectorWrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, 
		T & R,
		model_type & model, 
		time_type t)
  {
     // first eval RHS
    model.residual(*y.data(), *R.data(), t);
  }

private:
  friend residualPolicyBase<implicitEulerStandardResidual,
					 state_type,residual_type,
					 model_type, time_type, sizer_type>;

  friend advanceFullStatePolicyBase<implicitEulerStandardResidual,
				    state_type,residual_type,
				    model_type, time_type, sizer_type>;
  
};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
