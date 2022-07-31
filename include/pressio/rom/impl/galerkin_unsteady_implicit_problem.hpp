
#ifndef PRESSIO_ROM_IMPL_GALERKIN_IMPLICIT_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_IMPLICIT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class TrialSpaceType,
  class GalSystem,
  class StepperType
  >
class GalerkinImplicitProblem
{
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  GalSystem galSystem_;
  StepperType stepper_;

public:
  using independent_variable_type  = typename GalSystem::independent_variable_type;
  using state_type    = typename GalSystem::state_type;
  using residual_type = typename GalSystem::right_hand_side_type;
  using jacobian_type = typename GalSystem::jacobian_type;

  GalerkinImplicitProblem() = delete;

  template<class FomSystemType>
  GalerkinImplicitProblem(::pressio::ode::StepScheme schemeName,
			  TrialSpaceType & trialSpace,
			  const FomSystemType & fomSystem)
    : trialSpace_(trialSpace),
      galSystem_(trialSpace_, fomSystem),
      stepper_( ::pressio::ode::create_implicit_stepper(schemeName, galSystem_) )
  {}

  template<class ExtraArg>
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  ExtraArg && extra)
  {
    stepper_(state, sStart, sCount, sSize,
	     std::forward<ExtraArg>(extra));
  }

  residual_type createResidual() const{
    return stepper_.createResidual();
  }

  jacobian_type createJacobian() const{
    return stepper_.createJacobian();
  }

  void residual(const state_type & odeState, residual_type & R) const{
    stepper_.residual(odeState, R);
  }

  void jacobian(const state_type & odeState, jacobian_type & J) const{
    stepper_.jacobian(odeState, J);
  }
};

}}} // end pressio::rom::impl
#endif
