
#ifndef ODE_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"


namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type, 
	 typename sizer_type>
class implicitEulerJacobianStandardPolicy
  : public jacobianPolicyBase<implicitEulerJacobianStandardPolicy<
				state_type, jacobian_type,
				model_type, time_type, sizer_type> >
{
public:
  implicitEulerJacobianStandardPolicy() = default;
  ~implicitEulerJacobianStandardPolicy() = default;

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
		   time_type t,
		   time_type dt)
  {
    // first eval space jac
    model.jacobian( *y.data(), *J.data(), t);
    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
  }
  //----------------------------------------------------------------
  
private:
  friend jacobianPolicyBase<implicitEulerJacobianStandardPolicy<
			      state_type, jacobian_type,
			      model_type, time_type, sizer_type> >;

};//end class
  
}//end namespace polices
}//end namespace ode
#endif 
