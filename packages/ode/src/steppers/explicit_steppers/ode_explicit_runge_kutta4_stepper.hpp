
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_

#include "./impl/ode_explicit_runge_kutta4_stepper_impl.hpp"

namespace ode{

// namespace impl{
//   template <typename scalar>
//   struct btRK4 : public butcherTableau<scalar, 5, 5>{
//     // we put 5,5 so that we can store the coefficients with starting index 1
//     // to make it easier to match the theory
//     using base_t = butcherTableau<scalar, 5, 5>;
//     using base_t::a_;
//     using base_t::b_;
//     using base_t::c_;    
//     btRK4(){
//       a_(2,1) = static_cast<scalar>(1)/2;
//       a_(3,1) = static_cast<scalar>(0);
//       a_(4,1) = static_cast<scalar>(0);
//       a_(3,2) = static_cast<scalar>(1)/2;
//       a_(4,2) = static_cast<scalar>(0);
//       a_(4,3) = static_cast<scalar>(1);

//       c_[1] = static_cast<scalar>(0);
//       c_[2] = static_cast<scalar>(1)/2;
//       c_[3] = static_cast<scalar>(1)/2;
//       c_[4] = static_cast<scalar>(1);

//       b_[1] = static_cast<scalar>(1)/6;
//       b_[2] = static_cast<scalar>(1)/3;
//       b_[3] = static_cast<scalar>(1)/3;
//       b_[4] = static_cast<scalar>(1)/6;    
//     }
//   };
// }//end impl


  
///////////////////////
// Standard policy 
///////////////////////
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type
	 >
class ExplicitRungeKutta4Stepper<state_type,
				 residual_type,
				 scalar_type,
				 model_type,
				 time_type,
				 sizer_type,
				 void,
				 typename
				 std::enable_if<
				   !std::is_void<state_type>::value
				   >::type
				 >
  : public impl::ExplicitRungeKutta4StepperImpl<
	state_type, residual_type,
	scalar_type,
	model_type,
	time_type,
	sizer_type,
	ode::policy::explicit_residual_standard_policy<
	  state_type, residual_type,
	  model_type, time_type, sizer_type>>
{
public:
  using pol_t = ode::policy::explicit_residual_standard_policy<state_type,
						    residual_type,
						    model_type,
						    time_type,
						    sizer_type>;

  using base_t = impl::ExplicitRungeKutta4StepperImpl<state_type,
						      residual_type,
						      scalar_type,
						      model_type,
						      time_type,
						      sizer_type,
						      pol_t>;

public:
  template < typename T = model_type,
	     typename... Args>
  ExplicitRungeKutta4Stepper(T & model,
			     Args&&... rest)
    : base_t(model, policy_, std::forward<Args>(rest)...)
  {}

  ExplicitRungeKutta4Stepper() = delete;
  ~ExplicitRungeKutta4Stepper() = default;

private:
  pol_t policy_;

};//end class
  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
  
  
///////////////////////
// NON Standard policy 
///////////////////////
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type
	 >
class ExplicitRungeKutta4Stepper<state_type,
			   residual_type,
			   scalar_type,
			   model_type,
			   time_type,
			   sizer_type,
			   residual_policy_type,
			   typename
			   std::enable_if<
			     !std::is_void<residual_policy_type>::value
			     >::type
			   >
  : public impl::ExplicitRungeKutta4StepperImpl<state_type,
						residual_type,
						scalar_type,
						model_type,
						time_type,
						sizer_type,
						residual_policy_type>
{

public:
  using base_t = impl::ExplicitRungeKutta4StepperImpl<state_type,
						      residual_type,
						      scalar_type,
						      model_type,
						      time_type,
						      sizer_type,
						      residual_policy_type>;

public:
  template < typename T = model_type,
	     typename U = residual_policy_type,
	     typename... Args>
  ExplicitRungeKutta4Stepper(T & model,
		       U & policy,
		       Args&&... rest)
    : base_t(model, policy, std::forward<Args>(rest)...)
  {}
  ExplicitRungeKutta4Stepper() = delete;
  ~ExplicitRungeKutta4Stepper() = default;

};//end class

  
}//end namespace
#endif 








