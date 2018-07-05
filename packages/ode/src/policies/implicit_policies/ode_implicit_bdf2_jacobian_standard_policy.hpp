
#ifndef ODE_IMPLICIT_BDF2_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_BDF2_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "./impl/ode_bdf2_implicit_jacobian_impl.hpp"
#include "./base/ode_implicit_bdf2_jacobian_policy_base.hpp"
#include "../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type>
class implicitBDF2StandardJacobian
  : public implicitBDF2JacobianPolicyBase<implicitBDF2StandardJacobian,
					  state_type, jacobian_type,
					  model_type, time_type>,
    public advanceFullStatePolicyBase<implicitBDF2StandardJacobian,
				      state_type, jacobian_type,
				      model_type, time_type>
{
public:
  implicitBDF2StandardJacobian() = default;
  ~implicitBDF2StandardJacobian() = default;

private:
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
		   time_type t, time_type dt)
  {
    // first eval jac
    model.jacobian( *y.data(), *J.data(), t);

    jacobian_type A_( J.rows(),J.cols() );
    A_.setIdentity();

    // then fix it based on time stepping features
    ode::impl::implicit_euler_jacobian_impl(J, A_, dt);
  }

private:
  friend implicitBDF2JacobianPolicyBase<implicitBDF2StandardJacobian,
           state_type,jacobian_type,
           model_type, time_type>;

  friend advanceFullStatePolicyBase<implicitBDF2StandardJacobian,
				    state_type, jacobian_type,
				    model_type, time_type>;


};//end class
  
}//end namespace polices
}//end namespace ode
#endif 
