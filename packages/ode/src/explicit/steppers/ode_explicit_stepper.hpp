
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_traits.hpp"

namespace rompp{ namespace ode{

template<
  ExplicitEnum whichone,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename ...Args
  >
class ExplicitStepper
  : public ExplicitStepperBase<
  ExplicitStepper<
    whichone,
    state_type,
    system_type,
    velocity_type,
    Args...
    >
  >
{

  using this_t		= ExplicitStepper
    <whichone, state_type, system_type, velocity_type, Args...>;
  using base_t		= ExplicitStepperBase<this_t>;
  // need to friend base to allow it to access the () operator below
  friend base_t;

  using mytraits	= details::traits<this_t>;
  using scalar_type	= typename mytraits::scalar_t;
  using standard_res_policy_t = typename mytraits::standard_res_policy_t;
  using policy_t	= typename mytraits::velocity_policy_t;

  // this is the impl class type which holds all the implement details
  using impl_class_t	= typename mytraits::impl_t;
  impl_class_t myImpl_ = {};

  static constexpr auto zero = ::rompp::utils::constants::zero<scalar_type>();

public:
  ExplicitStepper()  = delete;
  ~ExplicitStepper() = default;

  // this is enabled all the time
  ExplicitStepper(state_type const	  & y0,
		  const system_type	  & model,
		  const policy_t	  & policyObj)
    : myImpl_(model,
	      policyObj,
	      y0,
	      policyObj(y0, model, this_t::zero)
	      )
  {}

  // only enable if the residual policy is standard
  template <
    typename T = policy_t,
    ::rompp::mpl::enable_if_t<
      mpl::is_same<
  	T, policy_t
  	>::value
      > * = nullptr
    >
  ExplicitStepper(const	state_type & y0,
  		  const system_type & model)
    : myImpl_(model,
  	      T(),
  	      y0,
  	      T()(y0, model, this_t::zero)
  	      )
  {}

private:
  // the compute method is private because we want users to use
  // the () operator in the base, which in turn calls compute here
  template<typename ... Args2>
  void compute(Args2 && ... args){
    myImpl_.doStep( std::forward<Args2>(args)... );
  }

};//end class

}} // end namespace rompp::ode
#endif
