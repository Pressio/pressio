
#ifndef ODE_IMPLICIT_ADAMS_MOULTON1_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_ADAMS_MOULTON1_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../impl/ode_adams_moulton1_implicit_residual_impl.hpp"
#include "../base/ode_implicit_adams_moulton1_residual_policy_base.hpp"
#include "../../common/ode_advance_full_state_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class implicitAdamsMoulton1StandardResidual
  : public implicitAdamsMoulton1ResidualPolicyBase<
  implicitAdamsMoulton1StandardResidual,
  state_type, residual_type, model_type, time_type>,
    public advanceFullStatePolicyBase<
  implicitAdamsMoulton1StandardResidual,
  state_type, residual_type, model_type, time_type>
{
public:
  implicitAdamsMoulton1StandardResidual() = default;
  ~implicitAdamsMoulton1StandardResidual() = default;  

private:

  void computeImpl(const state_type & y,
  		   const state_type & ynm1,
  		   const residual_type & RHSnm1,
  		   residual_type & R,
  		   model_type & model,
  		   time_type t,
  		   time_type dt)
  {
    if (R.empty())
      R.resize(y.size());
      
    // first eval RHS of the model
    this->computeModelResidual(y, R, model, t);
    // then fix residual based on time stepping features
    ode::impl::implicit_adams_moulton1_residual_impl(y, ynm1, RHSnm1, R, dt);
  }  
  //---------------------------------------------------
  
  void computeModelResidualImpl(const state_type & y,
				residual_type & RHS,
				model_type & model,
				time_type t)
  {
    if (RHS.empty())
      RHS.resize(y.size());

    model.residual(*y.data(), *RHS.data(), t);
  }  
  
private:
  friend implicitAdamsMoulton1ResidualPolicyBase<
  implicitAdamsMoulton1StandardResidual,
  state_type,residual_type,
  model_type, time_type>;

  friend advanceFullStatePolicyBase<
    implicitAdamsMoulton1StandardResidual,
    state_type,residual_type,
    model_type, time_type>;
  
};//end class
  
}//end namespace polices
}//end namespace ode  
#endif 
