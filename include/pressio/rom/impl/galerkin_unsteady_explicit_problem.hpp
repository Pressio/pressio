
#ifndef PRESSIO_ROM_IMPL_GALERKIN_EXPLICIT_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_EXPLICIT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <bool callableWithExtraArg, class GalSystem>
class GalerkinUnsteadyExplicitProblem
{
  using stepper_type =
    decltype(::pressio::ode::create_explicit_stepper(::pressio::ode::StepScheme::Euler,
						     std::declval<GalSystem &>()
						     ));
  GalSystem galSystem_;
  stepper_type stepper_;

public:
  using state_type = typename GalSystem::state_type;
  using independent_variable_type  = typename GalSystem::independent_variable_type;

  GalerkinUnsteadyExplicitProblem() = delete;

  template<class TrialSpaceType, class FomSystemType>
  GalerkinUnsteadyExplicitProblem(::pressio::ode::StepScheme schemeName,
				  const TrialSpaceType & trialSpace,
				  const FomSystemType & fomSystem)
    : galSystem_(trialSpace, fomSystem),
      stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_) )
  {}

  template<
    bool _callableWithExtraArg = callableWithExtraArg,
    mpl::enable_if_t<!_callableWithExtraArg, int> = 0>
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize)
  {
    stepper_(state, sStart, sCount, sSize);
  }

  template<
    class ExtraArg,
    bool _callableWithExtraArg = callableWithExtraArg,
    mpl::enable_if_t<_callableWithExtraArg, int> = 0
    >
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  ExtraArg && extra)
  {
    stepper_(state, sStart, sCount, sSize,
	     std::forward<ExtraArg>(extra));
  }
};

}}} // end pressio::rom::impl
#endif
