
#ifndef ODE_IMPLICIT_EULER_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../step_methods/impl/ode_euler_implicit_jacobian_impl.hpp"
#include "./base/ode_implicit_euler_jacobian_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type>
class implicitEulerStandardJacobian
  : public implicitEulerJacobianPolicyBase<
  implicitEulerStandardJacobian,state_type,
  jacobian_type,model_type,time_type>
{
public:
  implicitEulerStandardJacobian() = default;
  ~implicitEulerStandardJacobian() = default;
private:

  // // enable if using general type, not from core
  // template <typename U = state_type,
  // 	    typename T = jacobian_type,
  // 	    typename
  // 	    std::enable_if<
  // 	      core::meta::is_coreVectorWrapper<U>::value==false &&
  // 	      core::meta::is_coreMatrixWrapper<T>::value==false
  // 	      >::type * = nullptr
  // 	    >
  // void computeImpl(const U & y, T & J, model_type & model,
  // 		   time_type t, time_type dt){
  //   // first eval jac
  //   model.jacobian(y,J,t);
  //   // then fix it based on time stepping features
  //   ode::impl::implicit_euler_jacobian_impl(y, J, dt);
  // }
  // //---------------------------------------------------------

  // enable if using types from core package
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreMatrixWrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, T & J, model_type & model,
		   time_type t, time_type dt){
    // first eval jac
    model.jacobian( *y.data(), *J.data(), t);
    // then fix it based on time stepping features
    ode::impl::implicit_euler_jacobian_impl(J,dt);
  }
private:
  friend implicitEulerJacobianPolicyBase<implicitEulerStandardJacobian,
           state_type,jacobian_type,
           model_type, time_type>;

};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
