
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
    private OdeStorage<state_type, residual_type, 1>,
    private ImpOdeAuxData<model_type,
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
  using storage_base_t::auxStates_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;
  using auxdata_base_t::jacobian_obj_;
  using auxdata_base_t::t_;
  using auxdata_base_t::dt_;

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

    dt_ = dt;
    t_ = t;
    // store previous state = y;
    auxStates_[0] = y;
    y = solver.solve(*this, y);
  }//end doStepImpl
  //--------------------------------------------------------

protected:

  void residualImpl(const state_type & y, residual_type & R)const{
    (*residual_obj_).template operator()<my_enum, 1>(y, R, auxStates_, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  void jacobianImpl(const state_type & y, jacobian_type & J)const{
    (*jacobian_obj_).template operator()<my_enum>( y, J, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  residual_type residualImpl(const state_type & y)const{
    return (*residual_obj_).template operator()<my_enum, 1>( y, auxStates_, *model_, t_, dt_);
  }
  //--------------------------------------------------------

  jacobian_type jacobianImpl(const state_type & y)const{
    return (*jacobian_obj_).template operator()<my_enum>(y, *model_, t_, dt_);
  }
  //--------------------------------------------------------

private:
  friend stepper_base_t;

}; //end class

}}}//end namespace rompp::ode::impl
#endif
