
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_BDF2_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_BDF2_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename stateT,
	 typename residualT,
	 typename jacobianT,
	 typename model_type,
	 typename aux_stepper_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitBDF2StepperImpl<stateT,
			      residualT,
			      jacobianT,
			      model_type,
			      aux_stepper_type,
			      residual_policy_type,
			      jacobian_policy_type>
  : public ImplicitStepperBase<
             ImplicitBDF2StepperImpl<stateT,
				     residualT,
				     jacobianT,
				     model_type,
				     aux_stepper_type,
				     residual_policy_type,
				     jacobian_policy_type> >,
    protected OdeStorage<stateT, residualT, 2>,
    protected ImpOdeAuxData<model_type,
			  typename core::details::traits<stateT>::scalar_t,
			  residual_policy_type,
			  jacobian_policy_type>{

  static_assert( meta::is_legitimate_implicit_bdf2_residual_policy<
		 residual_policy_type>::value,
		 "IMPLICIT BDF2 RESIDUAL_POLICY NOT ADMISSIBLE,\
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  static_assert( meta::is_legitimate_implicit_bdf2_jacobian_policy<
		 jacobian_policy_type>::value,
		 "IMPLICIT BDF2 JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");


  using stepper_t = ImplicitEulerStepperImpl<stateT,
					     residualT,
					     jacobianT,
					     model_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using scalar_type  = typename core::details::traits<stateT>::scalar_t;
  using stepper_base_t = ImplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<stateT, residualT, 2>;
  using auxdata_base_t = ImpOdeAuxData<model_type, scalar_type,
				       residual_policy_type,
				       jacobian_policy_type>;
  static constexpr auto my_enum = ::rompp::ode::ImplicitEnum::BDF2;

public:
  using state_type = stateT;
  using residual_type = residualT;
  using jacobian_type = jacobianT;

public:
  template<typename solver_type, typename step_t>
  void operator()(stateT & y, scalar_type t,
		  scalar_type dt, step_t step,
		  solver_type & solver){
    this->dt_ = dt;
    this->t_ = t;

    // first step, use auxiliary stepper
    if (step == 1){
      this->auxStates_[0] = y;
      (*auxStp_)(y, t, dt, step, solver);
    }
    if (step == 2){
      this->auxStates_[1] = y;
      solver.solve(*this, y);
    }
    if (step >= 3){
      this->auxStates_[0] = this->auxStates_[1];
      this->auxStates_[1] = y;
      solver.solve(*this, y);
    }
  }//end doStepImpl

public:
  void residualImpl(const stateT & y, residualT & R)const{
    (*this->residual_obj_).template operator()<my_enum, 2>(y, R,
							   this->auxStates_,
							   *this->model_,
							   this->t_,
							   this->dt_);
  }
  //--------------------------------------------------------

  void jacobianImpl(const stateT & y, jacobianT & J)const{
    (*this->jacobian_obj_).template operator()<my_enum>(y, J,
							*this->model_,
							this->t_,
							this->dt_);
  }
  //--------------------------------------------------------

  residualT residualImpl(const stateT & y)const{
    return (*this->residual_obj_).template operator()<my_enum, 2>(y,
								  this->auxStates_,
								  *this->model_,
								  this->t_,
								  this->dt_);
  }
  //--------------------------------------------------------

  jacobianT jacobianImpl(const stateT & y)const{
    return (*this->jacobian_obj_).template operator()<my_enum>(y,
							       *this->model_,
							       this->t_,
							       this->dt_);
  }
  //--------------------------------------------------------

protected:
  aux_stepper_type * auxStp_;

protected:
  template < typename M = model_type,
       typename U = residual_policy_type,
       typename T = jacobian_policy_type,
       typename V = aux_stepper_type,
       typename T3 = stateT>
  ImplicitBDF2StepperImpl(const M & model,
        const U & res_policy_obj,
        const T & jac_policy_obj,
        V & auxStepper,
        const T3 & y0)
    : storage_base_t(y0),
      auxdata_base_t(model, res_policy_obj, jac_policy_obj),
      auxStp_(&auxStepper){}

  ImplicitBDF2StepperImpl() = delete;
  virtual ~ImplicitBDF2StepperImpl() = default;

private:
  friend stepper_base_t;


}; //end class

}}}//end namespace rompp::ode::impl
#endif
