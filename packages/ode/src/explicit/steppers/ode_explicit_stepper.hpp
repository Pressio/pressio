
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_traits.hpp"
#include "ode_explicit_stepper_helper_info.hpp"

namespace rompp{ namespace ode{

//!!!!!!!!!!!!!!!!
#ifdef HAVE_CPP14
//!!!!!!!!!!!!!!!!

    template<ExplicitEnum whichone, typename... Args>
    class ExplicitStepper
      : public impl::explicit_stepper_helper_info<whichone,
					 Args...>::base_impl_type{

      using info_t = impl::explicit_stepper_helper_info<whichone,
							Args...>;
      using state_t = typename info_t::state_type;
      using model_t = typename info_t::model_type;
      using res_t = typename info_t::res_type;
      using pol_t = typename info_t::residual_policy_type;
      using pol_std_t = typename info_t::res_std_pol_type;
      using base_impl_t = typename info_t::base_impl_type;

      //this needs to be public, it is detected by integrators
    public:
      using base_t = base_impl_t;

    public:

      // for standard policy, I need to have and pass the model
      ExplicitStepper(const model_t & model,
		      state_t const & y0,
		      res_t const & r0)
	: base_impl_t(model, policy_, y0, r0){}

      // for arbitrary policy, with also the model passed
      ExplicitStepper(const model_t & model,
		      const pol_t & policyObj,
		      state_t const & y0,
		      res_t const & r0)
	: base_impl_t(model, policyObj, y0, r0){}

      ExplicitStepper() = delete;
      ~ExplicitStepper() = default;

    private:
      // not used if policy is passed from outside
      pol_std_t policy_;
    };//end class


//!!!!!!!!!!!!!!!!!
#else
//!!!!!!!!!!!!!!!!!


    //if we have c++11 AND policy is STANDARD
    template<ExplicitEnum whichone,
	     typename ode_state_type,
	     typename model_type,
	     typename ode_residual_type>
    class ExplicitStepper<whichone, ode_state_type,
			  model_type, ode_residual_type, void>
      : public ExplicitStepperBase<
			ExplicitStepper<whichone, ode_state_type, model_type,
					ode_residual_type, void> >{

      using this_t	 = ExplicitStepper<whichone, ode_state_type,
					   model_type, ode_residual_type, void>;
      using mytraits	 = details::traits<this_t>;
      using base_t	 = ExplicitStepperBase<this_t>;
      friend base_t;
      using scalar_type  = typename mytraits::scalar_t;
      using pol_t	 = typename mytraits::residual_policy_t;
      using impl_class_t = typename details::traits<this_t>::impl_t;

      impl_class_t myImpl_ = {};

     public:
      ExplicitStepper(const model_type & model,
		      ode_state_type const & y0,
		      ode_residual_type const & r0)
	: myImpl_(model, pol_t(), y0, r0){}

      ExplicitStepper() = delete;
      ~ExplicitStepper() = default;

      template<typename step_t>
      void operator()(ode_state_type & y, scalar_type t,
		      scalar_type dt, step_t step){
	myImpl_.doStep(y, t, dt, step);
      }

    };//end class


    //if we have c++11 AND policy is user-defined
    template<ExplicitEnum whichone,
	     typename ode_state_type,
	     typename model_type,
	     typename ode_residual_type,
	     typename residual_policy_type>
    class ExplicitStepper<whichone, ode_state_type,
			  model_type, ode_residual_type,
			  residual_policy_type>
      : public ExplicitStepperBase<
		ExplicitStepper<whichone, ode_state_type, model_type,
				ode_residual_type, residual_policy_type> >{

      using this_t	 = ExplicitStepper<whichone, ode_state_type,
					   model_type, ode_residual_type,
					   residual_policy_type>;
      using mytraits	 = details::traits<this_t>;
      using base_t	 = ExplicitStepperBase<this_t>;
      friend base_t;
      using scalar_type  = typename mytraits::scalar_t;
      using impl_class_t = typename details::traits<this_t>::impl_t;

      impl_class_t myImpl_ = {};

     public:
      ExplicitStepper(const model_type & model,
		      const residual_policy_type & policyObj,
		      ode_state_type const & y0,
		      ode_residual_type const & r0)
	: myImpl_(model, policyObj, y0, r0){}

      ExplicitStepper() = delete;
      ~ExplicitStepper() = default;

      template<typename step_t>
      void operator()(ode_state_type & y, scalar_type t,
		      scalar_type dt, step_t step){
	myImpl_.doStep(y, t, dt, step);
      }
    };//end class

#endif

}} // end namespace rompp::ode
#endif
