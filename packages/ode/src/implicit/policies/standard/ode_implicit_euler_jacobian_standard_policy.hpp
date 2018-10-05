
#ifndef ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../ode_forward_declarations.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"

namespace rompp{ namespace ode{ namespace policy{

  
template<typename state_type,
	 typename jacobian_type,
	 typename model_type>
class ImplicitEulerJacobianStandardPolicy<
  state_type, jacobian_type, model_type,
  core::meta::enable_if_t<
    core::meta::is_core_vector_wrapper<state_type>::value and 
    core::meta::is_core_matrix_wrapper<jacobian_type>::value and
    std::is_same<typename core::details::traits<state_type>::scalar_t,
		 typename core::details::traits<jacobian_type>::scalar_t
		 >::value
    >
  >
  : public JacobianPolicyBase<ImplicitEulerJacobianStandardPolicy<
    state_type, jacobian_type,model_type> >{

public:
  ImplicitEulerJacobianStandardPolicy() = default;
  ~ImplicitEulerJacobianStandardPolicy() = default;

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;
    
private:
  
  // jacobian is passed by reference
  void operator()(const state_type & y, 
		  jacobian_type & J, 
		  model_type & model,
		  scalar_type t,
		  scalar_type dt){
    
    // first eval space jac
    model.jacobian( *y.data(), *J.data(), t);
    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(J, dt);
  }
  //----------------------------------------------------------------

  // jacobian is returned by the method
  jacobian_type operator()(const state_type & y, 
			   model_type & model,
			   scalar_type t,
			   scalar_type dt){
    
    auto nJJ = model.jacobian( *y.data(), t);
    jacobian_type JJ(nJJ);
    // update from time discrete residual
    ode::impl::implicit_euler_time_discrete_jacobian(JJ, dt);
    return JJ;
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
