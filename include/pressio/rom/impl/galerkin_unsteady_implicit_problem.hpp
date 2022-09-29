
#ifndef PRESSIO_ROM_IMPL_GALERKIN_IMPLICIT_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_IMPLICIT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class GalSystem>
class GalerkinUnsteadyImplicitProblem
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an explicit one
  using stepper_type =
    decltype(::pressio::ode::create_implicit_stepper
	     (::pressio::ode::StepScheme::BDF1,
	      std::declval<GalSystem &>()
	      ));

  using default_types = ImplicitGalerkinDefaultOperatorsTraits<typename GalSystem::state_type>;
  static_assert(std::is_same<
		typename default_types::reduced_residual_type,
		typename stepper_type::residual_type>::value,
		"galerkin implicit problem: impl: mismatching residual type");
  static_assert(std::is_same<
		typename default_types::reduced_jacobian_type,
		typename stepper_type::jacobian_type>::value,
		"galerkin implicit problem: impl: mismatching jacobian type");

public:
  using independent_variable_type  = typename GalSystem::independent_variable_type;
  using state_type = typename GalSystem::state_type;
  using residual_type = typename ImplicitGalerkinDefaultOperatorsTraits<state_type>::reduced_residual_type;
  using jacobian_type = typename ImplicitGalerkinDefaultOperatorsTraits<state_type>::reduced_jacobian_type;

  template<class ...Args>
  GalerkinUnsteadyImplicitProblem(::pressio::ode::StepScheme schemeName,
				  Args && ... args)
    : galSystem_(std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_implicit_stepper(schemeName, galSystem_) )
  {}

  stepper_type & galerkinStepper(){ return stepper_; }

  // template<class ExtraArg>
  // void operator()(state_type & state,
  // 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
  // 		  pressio::ode::StepCount sCount,
  // 		  pressio::ode::StepSize<independent_variable_type> sSize,
  // 		  ExtraArg && extra)
  // {
  //   stepper_(state, sStart, sCount, sSize,
  // 	     std::forward<ExtraArg>(extra));
  // }

private:
  GalSystem galSystem_;
  stepper_type stepper_;
};

}}} // end pressio::rom::impl
#endif
