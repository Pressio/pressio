
#ifndef ODE_CONCEPTS_OTHERS_HPP_
#define ODE_CONCEPTS_OTHERS_HPP_

namespace pressio{ namespace ode{

template <class T>
concept Steppable =
  requires(){
    typename T::independent_variable_type;
    typename T::state_type;
  }
  && requires(T & A,
	      typename T::state_type & state,
	      const ::pressio::ode::StepStartAt<typename T::independent_variable_type> & startAt,
	      const ::pressio::ode::StepCount & stepNumber,
	      const ::pressio::ode::StepSize<typename T::independent_variable_type> & dt)
  {
    A(state, startAt, stepNumber, dt);
  };

template <class T, class AuxT, class ...Args>
concept SteppableWithAuxiliaryArgs =
  requires(){
    typename T::independent_variable_type;
    typename T::state_type;
  }
  && requires(T & A,
	      typename T::state_type & state,
	      const ::pressio::ode::StepStartAt<typename T::independent_variable_type> & startAt,
	      const ::pressio::ode::StepCount & stepNumber,
	      const ::pressio::ode::StepSize<typename T::independent_variable_type> & dt,
	      AuxT && aux,
	      Args && ... args)
  {
    A(state, startAt, stepNumber, dt,
      std::forward<AuxT>(aux), std::forward<Args>(args)...);
  };

template <class T>
concept StronglySteppable = Steppable<T>;

template <class T, class AuxT, class ...Args>
concept StronglySteppableWithAuxiliaryArgs = SteppableWithAuxiliaryArgs<T, AuxT, Args...>;

template <class T, class IndVarType, class StateType>
concept StateObserver =
  requires(T && A,
	   const ::pressio::ode::StepCount & stepNumber,
	   const IndVarType & indVal,
	   const StateType & state)
  {
    A(stepNumber, indVal, state);
  };

template <class T, class IndVarType, class StateType>
concept StateGuesser =
  requires(T && A,
	   const ::pressio::ode::StepCount & stepNumber,
	   const ::pressio::ode::StepStartAt<IndVarType> & startAt,
	   StateType & state)
  {
    A(stepNumber, startAt, state);
  };

template <class T, class IndVarType>
concept StepSizePolicy =
  requires(T && A,
	   const ::pressio::ode::StepCount & stepNumber,
	   const ::pressio::ode::StepStartAt<IndVarType> & startAt,
	   ::pressio::ode::StepSize<IndVarType> & dt)
  {
    A(stepNumber, startAt, dt);
  };

template <class T, class IndVarType>
concept StepSizePolicyWithReductionScheme =
  requires(T && A,
	   const ::pressio::ode::StepCount & stepNumber,
	   const ::pressio::ode::StepStartAt<IndVarType> & startAt,
	   ::pressio::ode::StepSize<IndVarType> & dt,
	   ::pressio::ode::StepSizeMinAllowedValue<IndVarType> & minDt,
	   ::pressio::ode::StepSizeScalingFactor<IndVarType> & scalingFactor)
  {
    A(stepNumber, startAt, dt, minDt, scalingFactor);
  };

}}// end namespace pressio::ode
#endif
