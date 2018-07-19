
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_explicit_residual_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type>
class explicit_residual_standard_policy
  : public ExplicitResidualPolicyBase<
  explicit_residual_standard_policy<
    state_type, residual_type,
    model_type, time_type,
    sizer_type> >
{

public:
  explicit_residual_standard_policy() = default;
  ~explicit_residual_standard_policy() = default;  

private:
  //-----------------------------------------------------
  // enable when using types from core package
  //-----------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      core::meta::is_coreVector<U>::value==true &&
	      core::meta::is_coreVector<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y,
		   T & R,
		   model_type & model,
		   time_type t)
  {
    if (R.empty())
      sizer_type::matchSize(y, R);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
  }
  //-----------------------------------------------------

private:
  friend ExplicitResidualPolicyBase<
  explicit_residual_standard_policy<
    state_type, residual_type,
    model_type, time_type,
    sizer_type>>;

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
