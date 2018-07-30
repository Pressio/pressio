
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_explicit_residual_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename scalar_type,
	 typename sizer_type>
class explicit_residual_standard_policy
  : public ExplicitResidualPolicyBase<
  explicit_residual_standard_policy<
    state_type, residual_type,
    model_type, scalar_type,
    sizer_type> >
{

public:
  explicit_residual_standard_policy() = default;
  ~explicit_residual_standard_policy() = default;  

private:
  //-----------------------------------------------------
  // enable when state and residual are vector wrappers
  // what about the case when they are multivector wrappers?
  // think if this works right away
  //-----------------------------------------------------
  template <typename U = state_type,
	    typename T = residual_type,
	    typename
	    std::enable_if<
	      core::meta::is_core_vector<U>::value==true &&
	      core::meta::is_core_vector<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y,
		   T & R,
		   model_type & model,
		   scalar_type t)
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
    model_type, scalar_type,
    sizer_type>>;

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
