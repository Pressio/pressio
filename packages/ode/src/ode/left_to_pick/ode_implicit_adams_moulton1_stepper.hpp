
#ifndef ODE_IMPLICIT_ADAMS_MOULTON1_STEPPER_HPP_
#define ODE_IMPLICIT_ADAMS_MOULTON1_STEPPER_HPP_

#include "./impl/ode_implicit_adams_moulton1_stepper_impl.hpp"

namespace ode{

//*********************************************************
// residual policy = standard 
// jacobian policy = standard 
//*********************************************************
template<typename state_type, typename residual_type,
	 typename jacobian_type, typename scalar_type,
	 typename model_type, typename time_type,
	 typename sizer_type,
	 typename solver_policy_type
	 >
class implicitAdamsMoulton1Stepper<state_type, residual_type,
				   jacobian_type, scalar_type,
				   model_type, time_type, sizer_type,
				   solver_policy_type,
				   void, void,
				   typename std::enable_if<
				     !std::is_void<state_type>::value
				     >::type>
  : public impl::implicitAdamsMoulton1StepperImpl<
  state_type, residual_type, jacobian_type, scalar_type,
  model_type, time_type, sizer_type, solver_policy_type,
  ode::policy::residualStandardPolicy<
    state_type, residual_type,
    model_type, time_type, sizer_type>,
  ode::policy::jacobianStandardPolicy<
    state_type, jacobian_type,
    model_type, time_type, sizer_type>
  >
{

public:
  using res_pol_t = ode::policy::residualStandardPolicy<
  state_type, residual_type, model_type, time_type, sizer_type>;
			      
  using jac_pol_t = ode::policy::jacobianStandardPolicy<
  state_type, jacobian_type, model_type, time_type, sizer_type>;

  using base_t = impl::implicitAdamsMoulton1StepperImpl<state_type,
					       residual_type,
					       jacobian_type,
					       scalar_type,
					       model_type,
					       time_type,
					       sizer_type,
					       solver_policy_type,
					       res_pol_t,
					       jac_pol_t>;
public:
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename ... Args>
  implicitAdamsMoulton1Stepper(M & model,
		      S & solver,
		      Args&&... rest)
    : base_t(model, solver, res_policy_obj_, jac_policy_obj_,
    	std::forward<Args>(rest)...){}

  implicitAdamsMoulton1Stepper() = delete;
  ~implicitAdamsMoulton1Stepper() = default;

private:
  res_pol_t res_policy_obj_;
  jac_pol_t jac_policy_obj_;
  
}; //end class

  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
  

  
//*********************************************************
// residual policy = arbitrary, non void
// jacobian policy = arbitrary, non void
//*********************************************************
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type
	 >
class implicitAdamsMoulton1Stepper<state_type,
				   residual_type,
				   jacobian_type,
				   scalar_type,
				   model_type,
				   time_type,
				   sizer_type,
				   solver_policy_type,
				   residual_policy_type,
				   jacobian_policy_type,
				   typename
				   std::enable_if<
				     !std::is_void<residual_policy_type>::value &&
				     !std::is_void<jacobian_policy_type>::value
				     >::type
				   >
  : public impl::implicitAdamsMoulton1StepperImpl<state_type,
						  residual_type,
						  jacobian_type,
						  scalar_type,
						  model_type,
						  time_type,
						  sizer_type,
						  solver_policy_type,
						  residual_policy_type,
						  jacobian_policy_type>
{
public:
  using base_t = impl::implicitAdamsMoulton1StepperImpl<state_type,
							residual_type,
							jacobian_type,
							scalar_type,
							model_type,
							time_type,
							sizer_type,
							solver_policy_type,
							residual_policy_type,
							jacobian_policy_type>;
public:
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename ... Args>
  implicitAdamsMoulton1Stepper(M & model,
			       S & solver,
			       U & res_policy_obj,
			       T & jac_policy_obj,
			       Args&&... rest)
    : base_t(model, solver, res_policy_obj, jac_policy_obj,
	     std::forward<Args>(rest)...)
  {}
  implicitAdamsMoulton1Stepper() = delete;
  ~implicitAdamsMoulton1Stepper() = default;

}; //end class

  
}//end namespace ode
#endif
