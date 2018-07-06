
#ifndef ODE_EXPLICIT_RUNGE_KUTTA4_STANDARD_POLICY_HPP_
#define ODE_EXPLICIT_RUNGE_KUTTA4_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_explicit_residual_policy_base.hpp"
#include "../../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class explicitRungeKutta4StandardResidual
  : public explicitResidualPolicyBase<explicitRungeKutta4StandardResidual,
				      state_type, residual_type,
				      model_type, time_type,
				      sizer_type>,
    public advanceFullStatePolicyBase<explicitRungeKutta4StandardResidual,
				      state_type, residual_type,
				      model_type, time_type, sizer_type>
{
public:
  explicitRungeKutta4StandardResidual() = default;
  ~explicitRungeKutta4StandardResidual() = default;  

private:
  //----------------------------------------------------------------
  // enable when using types from core package
  //----------------------------------------------------------------
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
    if (R.empty())
      R.resize(y.size());

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
  }
  //----------------------------------------------------------------


  // void computeImpl(const state_type & y_prev,
  // 		   const state_type & y_curr,
  // 		   state_type & y_out,
  // 		   const residual_type & R1,
  // 		   const residual_type & R2,
  // 		   const residual_type & R3,
  // 		   const residual_type & R4,
  // 		   model_type & model,
  // 		   time_type t,
  // 		   time_type dt,
  // 		   time_type c1,
  // 		   time_type c2,
  // 		   time_type c3,
  // 		   time_type c4)
  // {
  //   model.residual(*y_curr.data(), *R4.data(), t);

  //   for (size_t i=0; i<y_inout.size(); i++){
  //     y_out[i] += c1*R1[i] + c2*R2[i] + c3*R3[i] + c4*R4[i];
  //   }    
  // }
  

private:
  friend explicitResidualPolicyBase<explicitRungeKutta4StandardResidual,
				    state_type, residual_type,
				    model_type, time_type, sizer_type>;

  friend advanceFullStatePolicyBase<explicitRungeKutta4StandardResidual,
				    state_type, residual_type,
				    model_type, time_type, sizer_type>;
  
};//end class

}//end namespace polices
}//end namespace ode  
#endif 
