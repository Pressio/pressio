
#ifndef ODE_EXPLICIT_EULER_STANDARD_POLICY_HPP_
#define ODE_EXPLICIT_EULER_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_explicit_policy_base.hpp"
#include "../../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class explicitEulerStandardResidual
  : public explicitResidualPolicyBase<explicitEulerStandardResidual,
				      state_type, residual_type,
				      model_type, time_type>,
    public advanceFullStatePolicyBase<explicitEulerStandardResidual,
				      state_type, residual_type,
				      model_type, time_type>

{
public:
  explicitEulerStandardResidual() = default;
  ~explicitEulerStandardResidual() = default;

private:
  //----------------------------------------------------------------
  // enable when using types from core package
  //----------------------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename sizer_t, 
	    typename
	    std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreVectorWrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, T & R,
		   model_type & model,
		   time_type t,
		   sizer_t & sizerObj)
  {
    if (R.empty())
      R.resize(y.size());

    model.residual(*y.data(), *R.data(), t);
  }
  //----------------------------------------------------------------

  // enable for general type, NOT from core
  template <typename U = state_type,
  	    typename T = residual_type,
	    typename sizer_t,
  	    typename
  	    std::enable_if<
  	      core::meta::is_coreVectorWrapper<U>::value==false &&
  	      core::meta::is_coreVectorWrapper<T>::value==false
  	      >::type * = nullptr
  	    >
  void computeImpl(const U & y, T & R,
		   model_type & model,
		   time_type t,
		   sizer_t & sizerObj)
  {
    if (sizerObj.getSize(R)==0)
      sizerObj.resize(R, sizerObj.getSize(y));

    model.residual(y, R, t);
  }

private:
  friend explicitResidualPolicyBase<explicitEulerStandardResidual,
				    state_type, residual_type,
				    model_type, time_type>;
  friend advanceFullStatePolicyBase<explicitEulerStandardResidual,
				    state_type, residual_type,
				    model_type, time_type>;


};//end class  

}//end namespace polices
}//end namespace ode  
#endif 
