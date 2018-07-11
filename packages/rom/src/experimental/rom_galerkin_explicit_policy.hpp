
#ifndef ROM_GALERKIN_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_RESIDUAL_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "policies/base/ode_advance_increment_policy_base.hpp"
#include "policies/base/ode_residual_policy_base.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class romGalerkinResidualPolicy
  : public ode::policy::residualPolicyBase<
  romGalerkinResidualPolicy, state_type, residual_type,
  model_type, time_type, sizer_type>,
    public ode::policy::advanceIncrementPolicyBase<
  romGalerkinResidualPolicy, state_type, residual_type,
  model_type, time_type, sizer_type>
{

private:
  using baseIncr_t = ode::policy::advanceIncrementPolicyBase<
  romGalerkinResidualPolicy, state_type, residual_type,
  model_type, time_type, sizer_type>;

public:
  romGalerkinResidualPolicy() = delete;
  romGalerkinResidualPolicy(const state_type & y0)
    : baseIncr_t(y0){}
  ~romGalerkinResidualPolicy() = default;  

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

private:
  // enable if using types from core package
  template <typename U = state_type, typename T = residual_type,
	    typename std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreVectorWrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y, 
		   T & R,
		   model_type & model, 
		   time_type t)
  {
    // reconstruct the solution
    yFull_ = *y0ptr_ + y;
    // eval RHS from target model
    R.setZero();
    model.residual(*yFull_.data(), *R.data(), t);
  }  

  // void weightTimeDiscreteResidualImpl(const state_type & y,
  // 				      residual_type & R,
  // 				      model_type & model, 
  // 				      time_type t)
  // {
  //   //no op
  // }
  
private:
  friend ode::policy::residualPolicyBase<
  romGalerkinResidualPolicy, state_type, residual_type,
  model_type, time_type, sizer_type>;

  friend baseIncr_t;

};//end class
  
}//end namespace exp
}//end namespace rom
#endif 
