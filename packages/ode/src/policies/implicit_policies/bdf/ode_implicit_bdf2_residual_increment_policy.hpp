
#ifndef ODE_IMPLICIT_BDF2_RESIDUAL_INCREMENT_POLICY_HPP_
#define ODE_IMPLICIT_BDF2_RESIDUAL_INCREMENT_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../impl/ode_bdf2_implicit_residual_impl.hpp"
#include "../base/ode_implicit_bdf2_residual_policy_base.hpp"
#include "../../common/ode_advance_increment_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class implicitBDF2IncrementResidual
  : public implicitBDF2ResidualPolicyBase<implicitBDF2IncrementResidual,
					   state_type, residual_type,
					   model_type, time_type>,
  public advanceIncrementPolicyBase<implicitBDF2IncrementResidual,
				    state_type, residual_type,
				    model_type, time_type>
{
private:
  using baseIncr_t = advanceIncrementPolicyBase<implicitBDF2IncrementResidual,
						state_type, residual_type,
						model_type, time_type>;
public:
  implicitBDF2IncrementResidual(const state_type & y0)
    : baseIncr_t(y0){}
  ~implicitBDF2IncrementResidual() = default;  

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

private:
  template <typename U = state_type, typename T = residual_type,
	    typename std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreVectorWrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y,
		   const U & ynm1,
		   const U & ynm2,
		   T & R,
		   model_type & model,
		   time_type t, time_type dt)
  { 
    // reconstruct the solution
    yFull_ = *y0ptr_ + y;

    // first eval RHS
    R.setZero();
    model.residual(*yFull_.data(), *R.data(), t);

    // then fix residual based on time stepping features
    ode::impl::implicit_bdf2_residual_impl(y, ynm1, ynm2, R, dt);
  }  
private:
  friend implicitBDF2ResidualPolicyBase<implicitBDF2IncrementResidual,
					 state_type, residual_type,
					 model_type, time_type>;
  friend baseIncr_t;

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
