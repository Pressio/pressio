
#ifndef ODE_INCREMENTBASED_JACOBIAN_HPP_
#define ODE_INCREMENTBASED_JACOBIAN_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../base/ode_advance_increment_policy_base.hpp"

namespace rompp{
namespace ode{
namespace policy{

template<typename state_type, typename jacobian_type,
	 typename model_type, typename time_type, typename sizer_type>
class incrementBasedJacobian
  : public ImplicitResidualPolicyBase<incrementBasedJacobian,
			      state_type, jacobian_type,
			      model_type, time_type, sizer_type>,
    public advanceIncrementPolicyBase<incrementBasedJacobian,
				      state_type, jacobian_type,
				      model_type, time_type, sizer_type>
{

private:
  using baseIncr_t = advanceIncrementPolicyBase<incrementBasedJacobian,
						state_type, jacobian_type,
						model_type, time_type, sizer_type>;

public:
  incrementBasedJacobian(const state_type & y0)
    : baseIncr_t(y0){}
  ~incrementBasedJacobian() = default;

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

private:
  // enable if using types from core package
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_core_matrix_wrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, T & J, model_type & model, time_type t)
  {
    // reconstruct the solution
    yFull_ = *y0ptr_ + y;
    // eval jac of the target model
    J.setZero();
    model.jacobian( *yFull_.data(), *J.data(), t);    
  }

  void weightJacobianImpl(const state_type & y,
			  jacobian_type & J,
			  model_type & model, 
			  time_type t)
  {
    // do nothing here
  }
  
private:
  friend ImplicitResidualPolicyBase<incrementBasedJacobian,
			    state_type,jacobian_type,
			    model_type, time_type, sizer_type>;

  friend baseIncr_t;

};//end class
  
}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 
