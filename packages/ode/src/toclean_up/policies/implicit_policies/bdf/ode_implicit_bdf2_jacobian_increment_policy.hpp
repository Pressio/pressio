
#ifndef ODE_IMPLICIT_BDF2_JACOBIAN_INCREMENT_POLICY_HPP_
#define ODE_IMPLICIT_BDF2_JACOBIAN_INCREMENT_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../impl/ode_bdf2_implicit_jacobian_impl.hpp"
#include "../base/ode_implicit_bdf2_jacobian_policy_base.hpp"
#include "../../common/ode_advance_increment_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename jacobian_type,
	 typename model_type, typename time_type>
class implicitBDF2IncrementJacobian
  : public implicitBDF2JacobianPolicyBase<implicitBDF2IncrementJacobian,
					   state_type, jacobian_type,
					   model_type, time_type>,
  public advanceIncrementPolicyBase<implicitBDF2IncrementJacobian,
			     state_type, jacobian_type,
			     model_type, time_type>
{
private:
  using baseIncr_t = advanceIncrementPolicyBase<implicitBDF2IncrementJacobian,
						state_type, jacobian_type,
						model_type, time_type>;

public:
  implicitBDF2IncrementJacobian(const state_type & y0)
    : baseIncr_t(y0){}
  ~implicitBDF2IncrementJacobian() = default;

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

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
    // reconstruct the solution
    yFull_ = *y0ptr_ + y;

    // first eval jac of the target model
    J.setZero();
    model.jacobian( *yFull_.data(), *J.data(), t);
    
    // fix jac based on euler backward
    jacobian_type A_( J.rows(), J.cols() );
    A_.setIdentity();
    // then fix it based on time stepping features
    ode::impl::implicit_bdf2_jacobian_impl(J, A_, dt);
  }
  
private:
  friend implicitBDF2JacobianPolicyBase<implicitBDF2IncrementJacobian,
           state_type,jacobian_type,
           model_type, time_type>;
  friend baseIncr_t;

};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
