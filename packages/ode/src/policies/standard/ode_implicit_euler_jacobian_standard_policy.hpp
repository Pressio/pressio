
#ifndef ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"


namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename sizer_type>
class implicit_euler_jacobian_standard_policy
  : public JacobianPolicyBase<implicit_euler_jacobian_standard_policy<
				state_type, jacobian_type,
				model_type, sizer_type> >
{
public:
  implicit_euler_jacobian_standard_policy() = default;
  ~implicit_euler_jacobian_standard_policy() = default;

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;
//   jacobian_type II_;
  
private:
  //----------------------------------------------------------------
  // enable if using types from core package
  //----------------------------------------------------------------
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename
	    std::enable_if<
	      core::meta::is_core_vector<U>::value==true &&
	      core::meta::is_coreMatrixWrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, 
		   T & J, 
		   model_type & model,
		   scalar_type t,
		   scalar_type dt)
  {
    // if (II_.rows()==0){
    //   II_.resize(J.rows(), J.cols());
    //   II_.setIdentity();
    // }
    // first eval space jac
    model.jacobian( *y.data(), *J.data(), t);
    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
  }
  //----------------------------------------------------------------
  
private:
  friend JacobianPolicyBase<
  implicit_euler_jacobian_standard_policy< state_type,
					   jacobian_type,
					   model_type,
					   sizer_type> >;

};//end class
  
}//end namespace polices
}//end namespace ode
#endif 
