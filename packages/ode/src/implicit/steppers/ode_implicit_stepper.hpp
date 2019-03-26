
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_

#include "ode_implicit_stepper_traits.hpp"
#include "ode_implicit_stepper_base.hpp"

namespace rompp{ namespace ode{

// both policies are STANDARD
template<ImplicitEnum whichone,
	 typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename aux_stepper_type>
class ImplicitStepper<whichone, ode_state_type,
		      ode_residual_type,
		      ode_jacobian_type, model_type,
		      aux_stepper_type, void, void>
  : public ImplicitStepperBase<
  ImplicitStepper<whichone, ode_state_type, ode_residual_type,
		  ode_jacobian_type, model_type,
		  aux_stepper_type, void, void>
  >{

  using this_t = ImplicitStepper<whichone, ode_state_type, ode_residual_type,
				 ode_jacobian_type, model_type,
				 aux_stepper_type, void, void>;
  using stepper_base_t = ImplicitStepperBase<this_t>;
  friend stepper_base_t;

  using mytraits	   = details::traits<this_t>;
  using residual_pol_t = typename mytraits::residual_policy_t;
  using jacobian_pol_t = typename mytraits::jacobian_policy_t;
  using scalar_t	   = typename mytraits::scalar_t;

  // this is the impl class
  using impl_class_t = typename mytraits::impl_t;
  impl_class_t myImpl = {};

public:
  // these need to be here because are detected by solver
  using state_type = ode_state_type;
  using residual_type = ode_residual_type;
  using jacobian_type = ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  // passing: model, initial state, and policies.
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model)
    : myImpl(model, residual_pol_t(), jacobian_pol_t(), y0){}

  ImplicitStepper(const ode_state_type & y0,
		  const ode_residual_type & r0,
		  const model_type & model)
    : myImpl(model, residual_pol_t(), jacobian_pol_t(), y0, r0){}

  // passing: model, initial state, aux_stepper
  // policy is standard for residual and jacobian
  template <typename T = aux_stepper_type,
	    core::meta::enable_if_t<
	      not std::is_void<T>::value
	      > * = nullptr>
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model,
		  T & auxStObj)
    : myImpl(model, residual_pol_t(),
	     jacobian_pol_t(), auxStObj, y0){}

public:
  template<typename step_t, typename solver_type,
	   typename ... Args>
  void operator()(ode_state_type & y,
		  scalar_t t,
		  scalar_t dt,
		  step_t step,
		  solver_type & solver,
		  Args ... args){
#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("current ode state is:","\n");
    ::rompp::core::io::print_core_wrapper(y, std::cout, 'd', -1);
#endif
    myImpl.doStep(*this, y, t, dt, step,
		  solver, std::forward<Args>(args)...);
  }

private:
  void residualImpl(const ode_state_type & y,
		    ode_residual_type & R) const{
    myImpl.computeResidual(y, R);
  }

  void jacobianImpl(const ode_state_type & y,
		    ode_jacobian_type & J) const{
    myImpl.computeJacobian(y, J);
  }

  ode_residual_type residualImpl(const ode_state_type & y) const{
    return 	myImpl.computeResidual(y);
  }

  ode_jacobian_type jacobianImpl(const ode_state_type & y) const{
    return myImpl.computeJacobian(y);
  }

};//end class




// residual is user-defined
// jacobian is standard
template<ImplicitEnum whichone,
	 typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename aux_stepper_type,
	 typename residual_policy_type>
class ImplicitStepper<whichone, ode_state_type, ode_residual_type,
		      ode_jacobian_type, model_type, aux_stepper_type,
		      residual_policy_type, void>
  : public ImplicitStepperBase<
  ImplicitStepper<whichone, ode_state_type, ode_residual_type,
		  ode_jacobian_type, model_type,
		  aux_stepper_type, residual_policy_type, void>
  >{

  using this_t = ImplicitStepper<whichone, ode_state_type, ode_residual_type,
				 ode_jacobian_type, model_type,
				 aux_stepper_type, residual_policy_type, void>;
  using stepper_base_t = ImplicitStepperBase<this_t>;
  friend stepper_base_t;

  using mytraits	   = details::traits<this_t>;
  using scalar_t	   = typename mytraits::scalar_t;
  using jacobian_pol_t = typename mytraits::jacobian_policy_t;

  // this is the impl class
  using impl_class_t = typename mytraits::impl_t;
  impl_class_t myImpl = {};

public:
  // these need to be here because are detected by solver
  using state_type = ode_state_type;
  using residual_type = ode_residual_type;
  using jacobian_type = ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  // passing: model, initial state, and policies
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj)
    : myImpl(model, resPolicyObj, jacobian_pol_t(), y0){}

  // passing: model, initial state,
  // initial residual and policies
  ImplicitStepper(const ode_state_type & y0,
		  const ode_residual_type & r0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj)
    : myImpl(model, resPolicyObj,
	     jacobian_pol_t(), y0, r0){}

  // passing: model, initial state, aux_stepper
  template <typename T = aux_stepper_type,
	    core::meta::enable_if_t<
	      not std::is_void<T>::value
	      > * = nullptr>
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj,
		  T & auxStObj)
    : myImpl(model, resPolicyObj,
	     jacobian_pol_t(),
	     auxStObj, y0){}

public:
  template<typename step_t, typename solver_type,
	   typename ... Args>
  void operator()(ode_state_type & y,
		  scalar_t t,
		  scalar_t dt,
		  step_t step,
		  solver_type & solver,
		  Args ... args){
#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("\ncurrent ode state is:","\n");
    ::rompp::core::io::print_core_wrapper(y, std::cout, 'd',-1);
#endif
    myImpl.doStep(*this, y, t, dt, step,
		  solver, std::forward<Args>(args)...);
  }

private:
  void residualImpl(const ode_state_type & y,
		    ode_residual_type & R) const{
    myImpl.computeResidual(y, R);
  }

  void jacobianImpl(const ode_state_type & y,
		    ode_jacobian_type & J) const{
    myImpl.computeJacobian(y, J);
  }

  ode_residual_type residualImpl(const ode_state_type & y) const{
    return 	myImpl.computeResidual(y);
  }

  ode_jacobian_type jacobianImpl(const ode_state_type & y) const{
    return myImpl.computeJacobian(y);
  }

};//end class







// both policies are user-defined
template<ImplicitEnum whichone,
	 typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename aux_stepper_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitStepper<whichone, ode_state_type, ode_residual_type,
		      ode_jacobian_type, model_type, aux_stepper_type,
		      residual_policy_type, jacobian_policy_type>
  : public ImplicitStepperBase<
  ImplicitStepper<whichone, ode_state_type, ode_residual_type,
		  ode_jacobian_type, model_type,
		  aux_stepper_type, residual_policy_type,
		  jacobian_policy_type>
  >{

  using this_t = ImplicitStepper<whichone, ode_state_type, ode_residual_type,
				 ode_jacobian_type, model_type,
				 aux_stepper_type, residual_policy_type,
				 jacobian_policy_type>;
  using stepper_base_t = ImplicitStepperBase<this_t>;
  friend stepper_base_t;

  using mytraits	   = details::traits<this_t>;
  using scalar_t	   = typename mytraits::scalar_t;

  // this is the impl class
  using impl_class_t = typename mytraits::impl_t;
  impl_class_t myImpl = {};

public:
  // these need to be here because are detected by solver
  using state_type = ode_state_type;
  using residual_type = ode_residual_type;
  using jacobian_type = ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  // passing: model, initial state, and policies
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj,
		  const jacobian_policy_type & jacPolicyObj)
    : myImpl(model, resPolicyObj, jacPolicyObj, y0){}

  // passing: model, initial state,
  // initial residual and policies
  ImplicitStepper(const ode_state_type & y0,
		  const ode_residual_type & r0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj,
		  const jacobian_policy_type & jacPolicyObj)
    : myImpl(model, resPolicyObj,
	     jacPolicyObj, y0, r0){}

  // passing: model, initial state, aux_stepper
  // policy is standard for residual and jacobian
  template <typename T = aux_stepper_type,
	    core::meta::enable_if_t<
	      not std::is_void<T>::value
	      > * = nullptr>
  ImplicitStepper(const ode_state_type & y0,
		  const model_type & model,
		  const residual_policy_type & resPolicyObj,
		  const jacobian_policy_type & jacPolicyObj,
		  T & auxStObj)
    : myImpl(model, resPolicyObj, jacPolicyObj,
	     auxStObj, y0){}

public:
  template<typename step_t, typename solver_type,
	   typename ... Args>
  void operator()(ode_state_type & y,
		  scalar_t t,
		  scalar_t dt,
		  step_t step,
		  solver_type & solver,
		  Args ... args){
#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("\ncurrent ode state is:","\n");
    ::rompp::core::io::print_core_wrapper(y, std::cout, 'd',-1);
#endif
    myImpl.doStep(*this, y, t, dt, step,
		  solver, std::forward<Args>(args)...);
  }

private:
  void residualImpl(const ode_state_type & y,
		    ode_residual_type & R) const{
    myImpl.computeResidual(y, R);
  }

  void jacobianImpl(const ode_state_type & y,
		    ode_jacobian_type & J) const{
    myImpl.computeJacobian(y, J);
  }

  ode_residual_type residualImpl(const ode_state_type & y) const{
    return 	myImpl.computeResidual(y);
  }

  ode_jacobian_type jacobianImpl(const ode_state_type & y) const{
    return myImpl.computeJacobian(y);
  }

};//end class


}} // end namespace rompp::ode
#endif
