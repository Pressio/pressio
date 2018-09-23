
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_RUNGEKUTTA4_STEPPER_HPP_

#include "./impl/ode_explicit_runge_kutta4_stepper_impl.hpp"

namespace rompp{
namespace ode{

/////////////////////////////////////////
//  (a) State and resid are same type
//  (b) Standard policy 
/////////////////////////////////////////

template<typename state_type,
	 typename model_type>
class ExplicitRungeKutta4Stepper<state_type, state_type,
				 model_type, void,
				 typename
				 std::enable_if<
				   !std::is_void<state_type>::value
				   >::type
				 >
  : public impl::ExplicitRungeKutta4StepperImpl<
	state_type, state_type,
	typename core::details::traits<state_type>::scalar_t,
        model_type,
	ode::policy::ExplicitResidualStandardPolicy<
	  state_type, state_type, model_type>>
{

public:
  using pol_t = ode::policy::ExplicitResidualStandardPolicy<
     state_type, state_type, model_type>;

  using base_t = impl::ExplicitRungeKutta4StepperImpl<
    state_type, state_type,
    typename core::details::traits<state_type>::scalar_t,
    model_type, pol_t>;

public:
  template < typename T1 = model_type,
	     typename T2 = state_type,
	     typename... Args>
  ExplicitRungeKutta4Stepper(T1 & model,
			     T2 const & y0,
			     T2 const & r0,
			     Args&&... rest)
    : base_t(model, policy_, y0, r0, std::forward<Args>(rest)...)
  {}

  ExplicitRungeKutta4Stepper() = delete;
  ~ExplicitRungeKutta4Stepper() = default;

private:
  pol_t policy_;

};//end class
  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
  

/////////////////////////////////////////
//  (a) State and resid are same type
//  (b) NON Standard policy 
/////////////////////////////////////////

template<typename state_type,
	 typename model_type,
	 typename residual_policy_type
	 >
class ExplicitRungeKutta4Stepper<state_type,
				 state_type,
				 model_type,
				 residual_policy_type,
				 typename
				 std::enable_if<
				   !std::is_void<residual_policy_type>::value
				   >::type
				 >
  : public impl::ExplicitRungeKutta4StepperImpl<
	     state_type,
	     state_type,
	     typename core::details::traits<state_type>::scalar_t,
	     model_type,
	     residual_policy_type>
{

public:
  using base_t = impl::ExplicitRungeKutta4StepperImpl<
	   state_type, state_type,
	   typename core::details::traits<state_type>::scalar_t,
	   model_type, residual_policy_type>;

public:
  template < typename T1 = model_type,
	     typename T2 = residual_policy_type,
	     typename T3 = state_type,
	     typename... Args>
 ExplicitRungeKutta4Stepper(T1 & model,
			    T2 & policy,
			    T3 const & y0,
			    T3 const & r0,
			    Args&&... rest)
    : base_t(model, policy, y0, r0,
	     std::forward<Args>(rest)...){}

  ExplicitRungeKutta4Stepper() = delete;
  ~ExplicitRungeKutta4Stepper() = default;

};//end class
  
}//end namespace
}//end namespace rompp
#endif 






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
