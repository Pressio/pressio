
#ifndef PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <int n, class = void> struct DeducedStepperType;

template <int n>
struct DeducedStepperType<n, mpl::enable_if_t< (n==-1) > >
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an implicit one
  template<class T>
  using type = decltype(::pressio::ode::create_implicit_stepper
			(::pressio::ode::StepScheme::BDF1,
			 std::declval<T &>()
			 ));
};

template <int n>
struct DeducedStepperType<n, mpl::enable_if_t< (n>0) > >
{
  template<class T>
  using type = decltype(::pressio::ode::create_implicit_stepper<n>(std::declval<T &>()));
};

template <
  class TrialSpaceType,
  class ResJacPolicyOrFullyDiscreteSystemType,
  int NumStates = -1
  >
class LspgUnsteadyProblem
{

  using stepper_type = typename DeducedStepperType<NumStates>::template
    type<ResJacPolicyOrFullyDiscreteSystemType>;

  using fom_states_manager_type = LspgFomStatesManager<TrialSpaceType>;

public:
  using independent_variable_type  = typename ResJacPolicyOrFullyDiscreteSystemType::independent_variable_type;
  using state_type    = typename TrialSpaceType::reduced_state_type;
  using residual_type = typename ResJacPolicyOrFullyDiscreteSystemType::residual_type;
  using jacobian_type = typename ResJacPolicyOrFullyDiscreteSystemType::jacobian_type;

  LspgUnsteadyProblem() = delete;

  template<
    class FomSystemType,
    class ...Args,
    int _NumStates = NumStates,
    mpl::enable_if_t< _NumStates == -1, int > = 0
    >
  LspgUnsteadyProblem(::pressio::ode::StepScheme odeSchemeName,
		      const TrialSpaceType & trialSpace,
		      const FomSystemType & fomSystem,
		      Args && ... args)
    : fomStatesManager_(create_lspg_fom_states_manager(odeSchemeName, trialSpace)),
      rjPolicyOrFullyDiscreteSystem_(trialSpace, fomSystem, fomStatesManager_, std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_implicit_stepper(odeSchemeName, rjPolicyOrFullyDiscreteSystem_))
  {}

  template<
    class FomSystemType,
    class ...Args,
    int _NumStates = NumStates,
    mpl::enable_if_t< (_NumStates > 0), int > = 0
    >
  LspgUnsteadyProblem(const TrialSpaceType & trialSpace,
		      const FomSystemType & fomSystem,
		      Args && ... args)
    : fomStatesManager_(create_lspg_fom_states_manager<_NumStates>(trialSpace)),
      rjPolicyOrFullyDiscreteSystem_(trialSpace, fomSystem, fomStatesManager_, std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_implicit_stepper<_NumStates>(rjPolicyOrFullyDiscreteSystem_))
  {}

  template<class SolverType, class ...ArgsOp>
  void operator()(state_type & reducedState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {
    stepper_(reducedState, sStart, sCount, sSize,
     	     solver, std::forward<ArgsOp>(argsop)...);
  }

  auto createState() const{
    return stepper_.createState();
  }

  residual_type createResidual() const{
    return stepper_.createResidual();
  }

  jacobian_type createJacobian() const{
    return stepper_.createJacobian();
  }

  void residualAndJacobian(const state_type & odeState,
			   residual_type & R,
			   jacobian_type & J,
			   bool computeJacobian) const{
    stepper_.residualAndJacobian(odeState, R, J, computeJacobian);
  }

private:
  LspgFomStatesManager<TrialSpaceType> fomStatesManager_;
  ResJacPolicyOrFullyDiscreteSystemType rjPolicyOrFullyDiscreteSystem_;
  stepper_type stepper_;
};

}}} // end pressio::rom::impl
#endif
