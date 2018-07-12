
#ifndef ODE_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../base/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type, 
	 typename sizer_type>
class jacobianStandardPolicy
  : public jacobianPolicyBase<jacobianStandardPolicy,
			      state_type, jacobian_type,
			      model_type, time_type, sizer_type>,
  public advanceFullStatePolicyBase<jacobianStandardPolicy,
				    state_type, jacobian_type,
				    model_type, time_type, sizer_type>
{
public:
  jacobianStandardPolicy() = default;
  ~jacobianStandardPolicy() = default;

private:
  //----------------------------------------------------------------
  // enable if using types from core package
  //----------------------------------------------------------------
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_coreMatrixWrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, 
		   T & J, 
		   model_type & model,
		   time_type t)
  {
    // first eval jac
    model.jacobian( *y.data(), *J.data(), t);

    jacobian_type A_( J.rows(),J.cols() );
    A_.setIdentity();
  }
  //----------------------------------------------------------------
  
private:
  friend jacobianPolicyBase<jacobianStandardPolicy,
           state_type,jacobian_type,
           model_type, time_type, sizer_type>;

  friend advanceFullStatePolicyBase<jacobianStandardPolicy,
				    state_type, jacobian_type,
				    model_type, time_type, sizer_type>;


};//end class
  
}//end namespace polices
}//end namespace ode
#endif 
