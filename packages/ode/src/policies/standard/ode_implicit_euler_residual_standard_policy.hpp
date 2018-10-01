
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../ode_ConfigDefs.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace rompp{
namespace ode{
namespace policy{

  
template<typename state_type, typename model_type,
	 typename residual_type = state_type, typename enable = void>
class ImplicitEulerResidualStandardPolicy;


//----------------------------------------------------------------
//----------------------------------------------------------------
//   default when state_type = residual_type
//----------------------------------------------------------------
//----------------------------------------------------------------
template<typename state_type, typename model_type>
class ImplicitEulerResidualStandardPolicy<
  state_type, model_type, state_type,
  core::meta::enable_if_t<
    core::meta::is_core_vector_wrapper<state_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitEulerResidualStandardPolicy<state_type, model_type>, 1, 0 >{
  
public:
  ImplicitEulerResidualStandardPolicy() = default;
  ~ImplicitEulerResidualStandardPolicy() = default;  

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;
  
  void computeImpl(const state_type & y, state_type & R,
  		   const std::array<state_type, 1> & oldYs,
  		   model_type & model, scalar_type t, scalar_type dt){
    
    if (R.empty())
      R.matchLayoutWith(y);

    R.setZero();
    model.residual(*y.data(), *R.data(), t);
    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
  }
  //----------------------------------------------------------------

  state_type computeImpl(const state_type & y, 
			 const std::array<state_type, 1> & oldYs,
			 model_type & model,
			 scalar_type t, scalar_type dt){
    
    auto nR = model.residual(*y.data(), t);
    state_type R(nR);
    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(y, oldYs[0], R, dt);
    return R;
  }
  //----------------------------------------------------------------
  
private:
  friend ImplicitResidualPolicyBase<
    ImplicitEulerResidualStandardPolicy<state_type, model_type,
					state_type>, 1,0>;
};//end class

}//end namespace polices
}//end namespace ode  
}//end namespace rompp
#endif 
