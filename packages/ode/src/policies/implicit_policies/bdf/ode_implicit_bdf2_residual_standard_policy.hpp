
#ifndef ODE_IMPLICIT_BDF2_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_BDF2_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../impl/ode_bdf2_implicit_residual_impl.hpp"
#include "../base/ode_implicit_bdf2_residual_policy_base.hpp"
#include "../../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class implicitBDF2StandardResidual
  : public implicitBDF2ResidualPolicyBase<implicitBDF2StandardResidual,
					  state_type, residual_type,
					  model_type, time_type>,
    public advanceFullStatePolicyBase<implicitBDF2StandardResidual,
				      state_type, residual_type,
				      model_type, time_type>
{
public:
  implicitBDF2StandardResidual() = default;
  ~implicitBDF2StandardResidual() = default;  

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
		   const U & ynm1,
		   const U & ynm2,
		   T & R,
		   model_type & model,
		   time_type t,
		   time_type dt)
  {
     // first eval RHS
    model.residual(*y.data(), *R.data(), t);
    // then fix residual based on time stepping features
    ode::impl::implicit_bdf2_residual_impl(y, ynm1, ynm2, R, dt);
  }  

private:
  friend implicitBDF2ResidualPolicyBase<implicitBDF2StandardResidual,
					 state_type,residual_type,
					 model_type, time_type>;

  friend advanceFullStatePolicyBase<implicitBDF2StandardResidual,
				    state_type,residual_type,
				    model_type, time_type>;
  
};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
