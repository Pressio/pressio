
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_

#include "ode_implicit_stepper_helper_info.hpp"
#include "ode_implicit_stepper_traits.hpp"
#include "./base/ode_implicit_stepper_base.hpp"

namespace rompp{ namespace ode{

//!!!!!!!!!!!!!!!!!
#if not defined HAVE_CPP14 // if NOT c++14, then we have c++11
//!!!!!!!!!!!!!!!!!

    //if we have c++11 AND policy is STANDARD
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
      ImplicitStepper(const model_type & model,
		      const ode_state_type & y0)
	: myImpl(model, residual_pol_t(), jacobian_pol_t(), y0){}

      ImplicitStepper(const model_type & model,
		      const ode_state_type & y0,
		      const ode_residual_type & r0)
	: myImpl(model, residual_pol_t(), jacobian_pol_t(), y0, r0){}

      // passing: model, initial state, aux_stepper
      // policy is standard for residual and jacobian
      template <typename T = aux_stepper_type,
		core::meta::enable_if_t<
		  not std::is_void<T>::value
		  > * = nullptr>
      ImplicitStepper(const model_type & model,
		      const ode_state_type & y0,
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



    //if we have c++11 AND policy is user defined
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

#ifdef DEBUG_PRINT
      int statePrintSize_ = -1;
#endif

    public:
      ImplicitStepper() = delete;
      ~ImplicitStepper() = default;

#ifdef DEBUG_PRINT
      void setStatePrintSize(int n){
	statePrintSize_ = n;
      }
#endif

      // passing: model, initial state, and policies
      ImplicitStepper(const model_type & model,
    		      const residual_policy_type & resPolicyObj,
    		      const jacobian_policy_type & jacPolicyObj,
    		      const ode_state_type & y0)
    	: myImpl(model, resPolicyObj, jacPolicyObj, y0){}

      // passing: model, initial state,
      // initial residual and policies
      ImplicitStepper(const model_type & model,
    		      const residual_policy_type & resPolicyObj,
    		      const jacobian_policy_type & jacPolicyObj,
    		      const ode_state_type & y0,
    		      const ode_residual_type & r0)
    	: myImpl(model, resPolicyObj,
		 jacPolicyObj, y0, r0){}

      // passing: model, initial state, aux_stepper
      // policy is standard for residual and jacobian
      template <typename T = aux_stepper_type,
    		core::meta::enable_if_t<
    		  not std::is_void<T>::value
    		  > * = nullptr>
      ImplicitStepper(const model_type & model,
    		      const residual_policy_type & resPolicyObj,
    		      const jacobian_policy_type & jacPolicyObj,
    		      const ode_state_type & y0,
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
	::rompp::core::io::print_core_wrapper(y, std::cout, 'd',
					       statePrintSize_);
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



// //!!!!!!!!!!!!!!!!!
// #else // if we have C++14
// //!!!!!!!!!!!!!!!!!

// TODO: revise this for when c++14 is back

//     template<ImplicitEnum whichone, typename... Args>
//     class ImplicitStepper
//       : public impl::implicit_stepper_helper_info<whichone,
// 						  Args...>::base_impl_type{

//       using info_t = impl::implicit_stepper_helper_info<whichone,
// 							Args...>;
//       using state_t = typename info_t::state_type;
//       using model_t = typename info_t::model_type;
//       using res_t = typename info_t::res_type;
//       using jac_t = typename info_t::jac_type;
//       using residual_pol_t = typename info_t::residual_policy_type;
//       using residual_pol_std_t = typename info_t::res_std_pol_type;
//       using jacobian_pol_t = typename info_t::jacobian_policy_type;
//       using jacobian_pol_std_t = typename info_t::jac_std_pol_type;
//       using auxiliary_stepper_t = typename info_t::auxiliary_stepper_type;
//       using base_impl_t = typename info_t::base_impl_type;

//       //this needs to be public, it is detected by integrators
//     public:
//       using base_t = base_impl_t;

//     public:
//       // some constructors need y0, and some need y0/r0
//       // because some multi-stage methods need history of RHS
//       // e.g. backward Euler only needs y, y_n-1
//       // e.g. BDF2 only needs y, y_n-1, y_n-2

//       // passing: model, initial state
//       // policy is standard for residual and jacobian
//       ImplicitStepper(const model_t & model,
// 		      const state_t & y0)
// 	: base_impl_t(model, res_policy_obj_,
// 		      jac_policy_obj_, y0){}

//       // passing: model, initial state, aux_stepper
//       // policy is standard for residual and jacobian
//       template <typename T = auxiliary_stepper_t,
// 		core::meta::enable_if_t<
// 		  not std::is_void<T>::value
// 		  > * = nullptr>
//       ImplicitStepper(const model_t & model,
// 		      const state_t & y0,
// 		      T & auxStObj)
// 	: base_impl_t(model, res_policy_obj_,
// 		      jac_policy_obj_, auxStObj, y0){}

//       // passing: model, initial state, and policies.
//       ImplicitStepper(const model_t & model,
// 		      const residual_pol_t & resPolicyObj,
// 		      const jacobian_pol_t & jacPolicyObj,
// 		      const state_t & y0)
// 	: base_impl_t(model, resPolicyObj, jacPolicyObj, y0){}

//       // passing: model, initial state, initial residual
//       // policy is standard for residual and jacobian
//       ImplicitStepper(const model_t & model,
// 		      const state_t & y0,
// 		      const res_t & r0)
// 	: base_impl_t(model, res_policy_obj_,
// 		      jac_policy_obj_, y0, r0){}

//       // passing: model, initial state,
//       // initial residual and policies.
//       ImplicitStepper(const model_t & model,
// 		      const residual_pol_t & resPolicyObj,
// 		      const jacobian_pol_t & jacPolicyObj,
// 		      const state_t & y0,
// 		      const res_t & r0)
// 	: base_impl_t(model, resPolicyObj,
// 		      jacPolicyObj, y0, r0){}

//       ImplicitStepper() = delete;
//       virtual ~ImplicitStepper() = default;

//     private:
//       // not used if policy is passed from outside
//       residual_pol_std_t res_policy_obj_;
//       jacobian_pol_std_t jac_policy_obj_;
//     };//end class

#endif

}} // end namespace rompp::ode
#endif
