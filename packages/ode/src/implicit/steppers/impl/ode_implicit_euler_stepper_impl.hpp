
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename model_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitEulerStepperImpl<state_type,
			       residual_type,
			       jacobian_type,
			       model_type,
			       residual_policy_type,
			       jacobian_policy_type>
  : public ImplicitStepperBase<
	    ImplicitEulerStepperImpl<state_type, residual_type,
				     jacobian_type,
				     model_type,
				     residual_policy_type,
				     jacobian_policy_type> >,
    protected OdeStorage<state_type, residual_type, 1>,
    protected ImpOdeAuxData<model_type,
			  typename core::details::traits<state_type>::scalar_t,
			  residual_policy_type, jacobian_policy_type>{

  static_assert( meta::is_legitimate_implicit_euler_residual_policy<
		 residual_policy_type>::value,
"IMPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  static_assert( meta::is_legitimate_implicit_euler_jacobian_policy<
		 jacobian_policy_type>::value,
"IMPLICIT EULER JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  using stepper_t = ImplicitEulerStepperImpl<state_type,
					     residual_type,
					     jacobian_type,
					     model_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using scalar_type  = typename core::details::traits<state_type>::scalar_t;
  using stepper_base_t = ImplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, residual_type, 1>;
  using auxdata_base_t = ImpOdeAuxData<model_type, scalar_type,
				       residual_policy_type,
				       jacobian_policy_type>;

  static constexpr auto my_enum = ::rompp::ode::ImplicitEnum::Euler;

public:
  // these aliases are needed by the solver
  using vector_type = state_type;
  using matrix_type = jacobian_type;

protected:
  template < typename M = model_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename T3 = state_type>
  ImplicitEulerStepperImpl(const M & model,
			   const U & res_policy_obj,
			   const T & jac_policy_obj,
			   const T3 & y0)
    : storage_base_t(y0),
      auxdata_base_t(model, res_policy_obj, jac_policy_obj){}

  ImplicitEulerStepperImpl() = delete;
  virtual ~ImplicitEulerStepperImpl(){};

public:
  template<typename solver_type, typename step_t>
  void operator()(state_type & y, scalar_type t,
		  scalar_type dt, step_t step,
		  solver_type & solver){

    this->dt_ = dt;
    this->t_ = t;
    // store previous state = y;
    this->auxStates_[0] = y;
    solver.solve(*this, y);
  }//end doStepImpl
  //--------------------------------------------------------

public:
  void residualImpl(const state_type & y, residual_type & R)const{
    (*this->residual_obj_).template operator()<my_enum, 1>(y, R,
							   this->auxStates_,
							   *this->model_,
							   this->t_,
							   this->dt_);
  }
  //--------------------------------------------------------

  void jacobianImpl(const state_type & y, jacobian_type & J)const{
    (*this->jacobian_obj_).template operator()<my_enum>( y, J,
							 *this->model_,
							 this->t_,
							 this->dt_);
  }
  //--------------------------------------------------------

  residual_type residualImpl(const state_type & y)const{
    return (*this->residual_obj_).template operator()<my_enum, 1>( y,
								   this->auxStates_,
								   *this->model_,
								   this->t_,
								   this->dt_);
  }
  //--------------------------------------------------------

  jacobian_type jacobianImpl(const state_type & y)const{
    return (*this->jacobian_obj_).template operator()<my_enum>(y,
							       *this->model_,
							       this->t_,
							       this->dt_);
  }
  //--------------------------------------------------------

private:
  friend stepper_base_t;

}; //end class

}}}//end namespace rompp::ode::impl
#endif
