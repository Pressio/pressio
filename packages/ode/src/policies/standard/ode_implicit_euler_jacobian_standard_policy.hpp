
#ifndef ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../ode_ConfigDefs.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"

namespace rompp{
namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type>
class ImplicitEulerJacobianStandardPolicy
  : public JacobianPolicyBase<ImplicitEulerJacobianStandardPolicy<
				state_type, jacobian_type,
				model_type> >
{
public:
  ImplicitEulerJacobianStandardPolicy() = default;
  ~ImplicitEulerJacobianStandardPolicy() = default;

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;
  
private:
  template <typename U = state_type,
	    typename T = jacobian_type,
	    typename std::enable_if<
	      core::meta::is_core_vector_wrapper<U>::value==true &&
	      core::meta::is_core_matrix_wrapper<T>::value==true
	      >::type * = nullptr
	    >
  void computeImpl(const U & y, 
		   T & J, 
		   model_type & model,
		   scalar_type t,
		   scalar_type dt)
  {

    // first eval space jac
    model.jacobian( *y.data(), *J.data(), t);
    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
  }
  
private:
  friend JacobianPolicyBase<
  ImplicitEulerJacobianStandardPolicy< state_type,
					   jacobian_type,
					   model_type> >;

};//end class
  
}//end namespace polices
}//end namespace ode
}//end namespace rompp
#endif 
