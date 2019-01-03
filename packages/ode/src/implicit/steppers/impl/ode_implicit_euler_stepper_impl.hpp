
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename stateT,
	 typename residualT,
	 typename jacobianT,
	 typename model_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitEulerStepperImpl<stateT,
			       residualT,
			       jacobianT,
			       model_type,
			       residual_policy_type,
			       jacobian_policy_type>
  : protected OdeStorage<stateT, residualT, 1>,
    protected ImpOdeAuxData<model_type,
			  typename core::details::traits<stateT>::scalar_t,
			  residual_policy_type, jacobian_policy_type>
{

  static_assert( meta::is_legitimate_implicit_euler_residual_policy<
		 residual_policy_type>::value,
"IMPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  static_assert( meta::is_legitimate_implicit_euler_jacobian_policy<
		 jacobian_policy_type>::value,
"IMPLICIT EULER JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OR DERIVING FROM WRONG BASE");

  using this_t = ImplicitEulerStepperImpl<stateT, residualT,
					  jacobianT, model_type,
					  residual_policy_type,
					  jacobian_policy_type>;
  using scalar_type  = typename core::details::traits<stateT>::scalar_t;
  using storage_base_t = OdeStorage<stateT, residualT, 1>;
  using auxdata_base_t = ImpOdeAuxData<model_type, scalar_type,
				       residual_policy_type,
				       jacobian_policy_type>;

  static constexpr auto my_enum = ::rompp::ode::ImplicitEnum::Euler;

public:
  ImplicitEulerStepperImpl(const model_type & model,
			   const residual_policy_type & res_policy_obj,
			   const jacobian_policy_type & jac_policy_obj,
			   const stateT & y0)
    : storage_base_t(y0),
      auxdata_base_t(model, res_policy_obj, jac_policy_obj){}

  ImplicitEulerStepperImpl() = delete;
  ~ImplicitEulerStepperImpl(){};

public:
  template<typename solver_type, typename step_t,
	   typename caller_t>
  void doStep(stateT & y, scalar_type t,
		scalar_type dt, step_t step,
		solver_type & solver,
		caller_t & caller){

    this->dt_ = dt;
    this->t_ = t;
    this->auxStates_[0] = y;
    solver.solve(caller, y);
  }

  void computeResidual(const stateT & y, residualT & R) const{
    (*this->residual_obj_).template operator()<my_enum, 1>(y, R,
							   this->auxStates_,
							   *this->model_,
							   this->t_,
							   this->dt_);
  }

  void computeJacobian(const stateT & y, jacobianT & J) const{
    (*this->jacobian_obj_).template operator()<my_enum>( y, J,
							 *this->model_,
							 this->t_,
							 this->dt_);
  }

  residualT computeResidual(const stateT & y) const{
    return (*this->residual_obj_).template operator()<my_enum, 1>( y,
								   this->auxStates_,
								   *this->model_,
								   this->t_,
								   this->dt_);
  }

  jacobianT computeJacobian(const stateT & y) const{
    return (*this->jacobian_obj_).template operator()<my_enum>(y,
							       *this->model_,
							       this->t_,
							       this->dt_);
  }

// private:
//   friend stepper_base_t;

}; //end class

}}}//end namespace rompp::ode::impl
#endif
