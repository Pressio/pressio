
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../ode_ConfigDefs.hpp"
#include "../base/ode_explicit_residual_policy_base.hpp"

namespace rompp{
namespace ode{
namespace policy{

template<typename state_type,
	 typename space_residual_type,
	 typename model_type>
class ExplicitResidualStandardPolicy
  : public ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<
    state_type, space_residual_type, model_type> >
{

public:
  ExplicitResidualStandardPolicy() = default;
  ~ExplicitResidualStandardPolicy() = default;  
  
private:

  //-----------------------------------------------------
  // enable when state and residual are vector wrappers
  // what about the case when they are multivector wrappers?
  // think if this works right away
  //-----------------------------------------------------
  template <typename U = state_type,
	    typename T = space_residual_type,
	    typename sc_t,
	    typename
	    std::enable_if<
	      core::meta::is_core_vector_wrapper<U>::value==true &&
	      std::is_same<U, T>::value &&
	      std::is_same<typename core::details::traits<U>::scalar_t,
			   sc_t>::value
	      /*core::meta::is_core_multi_vector<U>::value==true &&
	      core::meta::is_core_multi_vector<T>::value==true*/
	      >::type * = nullptr
	    >
  void operator()(const U & y, T & R, model_type & model, sc_t t)
  {
    if (R.empty())
      R.matchLayoutWith(y);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
  }
  //-----------------------------------------------------

private:
  friend ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<
    state_type, space_residual_type, model_type> >;

};//end class

}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 
