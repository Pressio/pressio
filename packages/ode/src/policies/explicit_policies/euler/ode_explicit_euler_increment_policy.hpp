
#ifndef ODE_EXPLICIT_EULER_INCREMENT_POLICY_HPP_
#define ODE_EXPLICIT_EULER_INCREMENT_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_explicit_residual_policy_base.hpp"
#include "../../common/ode_advance_increment_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class explicitEulerIncrementResidual
  : public explicitResidualPolicyBase<explicitEulerIncrementResidual,
				      state_type, residual_type,
				      model_type, time_type, sizer_type>,
    public advanceIncrementPolicyBase<explicitEulerIncrementResidual,
				      state_type, residual_type,
				      model_type, time_type, sizer_type>
{
private:
  using baseIncr_t = advanceIncrementPolicyBase<explicitEulerIncrementResidual,
						state_type, residual_type,
						model_type, time_type>;
public:
  explicitEulerIncrementResidual(const state_type & y0) : baseIncr_t(y0){}
  ~explicitEulerIncrementResidual() = default;  

private:
  using baseIncr_t::yFull_;
  using baseIncr_t::y0ptr_;

private:
  //----------------------------------------------------------------
  // enable when using types from core package
  //----------------------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename std::enable_if<
	      core::meta::is_coreVectorWrapper<U>::value==true &&
	      core::meta::is_coreVectorWrapper<T>::value==true
	    >::type * = nullptr>
  void computeImpl(const U & y, T & R, model_type & model, time_type t)
  { 
    if (R.empty())
      R.resize(y.size());

    // reconstruct the solution
    yFull_ = *y0ptr_ + y;
    // eval RHS
    R.setZero();
    model.residual(*yFull_.data(), *R.data(), t);
  }  

  // //----------------------------------------------------------------
  // // enable for general type, NOT from core
  // //----------------------------------------------------------------
  // template <typename U = state_type,
  // 	    typename T = residual_type,
  // 	    typename
  // 	    std::enable_if<
  // 	      core::meta::is_coreVectorWrapper<U>::value==false &&
  // 	      core::meta::is_coreVectorWrapper<T>::value==false
  // 	      >::type * = nullptr
  // 	    >
  // void computeImpl(const U & y, T & R, model_type & model, time_type t)
  // {
  //   auto stateSz = sizer_type::getSize(y);
  //   if (sizer_type::getSize(R)==0)
  //     sizer_type::resize(R, sizer_type::getSize(y));

  //   // reconstruct the solution
  //   for (size_t i=0; i<stateSz; i++)
  //     yFull_[i] = *y0ptr_[i] + y[i];

  //   auto zero = static_cast<decltype(R[0])>(0);
  //   for (size_t i=0; i<stateSz; i++)
  //     R[i] = zero;
    
  //   model.residual(yFull_, R, t);
  // }

private:
  friend explicitResidualPolicyBase<explicitEulerIncrementResidual,
				    state_type, residual_type,
				    model_type, time_type,
				    sizer_type>;
  friend baseIncr_t;

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
