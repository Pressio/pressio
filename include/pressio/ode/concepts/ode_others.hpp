
#ifndef ODE_CONCEPTS_OTHERS_CXX17_HPP_
#define ODE_CONCEPTS_OTHERS_CXX17_HPP_

namespace pressio{ namespace ode{ namespace impl{

/*
  we need to ensure operator() accepts the state
  by **non-const** reference. One (maybe only) way to detect
  that is check that if operator() can bind an rvalue.
  Basically if we tried to bind an rvalue to an lvalue reference,
  it should fail, and if it fails that is good because it is what we want.
*/

template <class T, class IndVarType, class StateType, class Enable = void>
struct guesser_taking_state_by_ref
  : std::true_type{};

template <class T, class IndVarType, class StateType>
struct guesser_taking_state_by_ref<
  T, IndVarType, StateType,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< StateType >()
	)
       )
      >::value
    >
  > : std::false_type{};

/*
  we need to ensure operator() accepts the step size
  by **non-const** reference. One (maybe only) way to detect
  that is check that if operator() can bind an rvalue.
  Basically if we tried to bind an rvalue to an lvalue reference,
  it should fail, and if it fails that is good because it is what we want.
*/

template <class T, class IndVarType, class Enable = void>
struct step_size_policy_taking_dt_by_ref
  : std::true_type{};

template <class T, class IndVarType>
struct step_size_policy_taking_dt_by_ref<
  T, IndVarType,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> >()
	)
       )
      >::value
    >
  > : std::false_type{};


#define STEP_SIZE_POLICY_WITH_REDUC_TAKING_DT_BY_REF(T1, T2, T3)	\
  T, IndVarType, \
  mpl::enable_if_t<							\
    std::is_void<							\
    decltype								\
  (									\
 std::declval<T>()							\
   (									\
    std::declval< ::pressio::ode::StepCount >(),			\
      std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),	\
      std::declval< T1 >(),						\
      std::declval< T2 >(),						\
      std::declval< T3 >()						\
      )									\
   )									\
  >::value \
 > \

template <class T, class IndVarType, class Enable = void>
struct step_size_policy_with_reduc_taking_dt_by_ref
  : std::true_type{};

template <class T, class IndVarType>
struct step_size_policy_with_reduc_taking_dt_by_ref<
  STEP_SIZE_POLICY_WITH_REDUC_TAKING_DT_BY_REF(StepSize<IndVarType>, \
						    StepSizeMinAllowedValue<IndVarType>, \
						    StepSizeScalingFactor<IndVarType>)
  > : std::false_type{};

template <class T, class IndVarType>
struct step_size_policy_with_reduc_taking_dt_by_ref<
  STEP_SIZE_POLICY_WITH_REDUC_TAKING_DT_BY_REF(StepSize<IndVarType>, \
						    StepSizeMinAllowedValue<IndVarType> &, \
						    StepSizeScalingFactor<IndVarType> &)
  > : std::false_type{};

template <class T, class IndVarType>
struct step_size_policy_with_reduc_taking_dt_by_ref<
  STEP_SIZE_POLICY_WITH_REDUC_TAKING_DT_BY_REF(StepSize<IndVarType> &, \
						    StepSizeMinAllowedValue<IndVarType> &, \
						    StepSizeScalingFactor<IndVarType>)
  > : std::false_type{};

template <class T, class IndVarType>
struct step_size_policy_with_reduc_taking_dt_by_ref<
  STEP_SIZE_POLICY_WITH_REDUC_TAKING_DT_BY_REF(StepSize<IndVarType> &, \
						    StepSizeMinAllowedValue<IndVarType> , \
						    StepSizeScalingFactor<IndVarType> &)
  > : std::false_type{};

// we are missing all cases with just a single one that is ref because
// when i tried to add those special cases I get ambiguous template error

/*
  for any stepper we need to ensure its operator() accepts the
  state by **non-const** reference. One (maybe only) way to detect
  that is check that if operator() can bind an rvalue state.
  Basically if we tried to bind an rvalue to an lvalue reference,
  it should fail, and if it fails that is good because it is what we want.

  In brief, we want the stepper to have this:
  class{
      void operator()(state_type & state, ...)
  };

  NOT this:
  class{
      void operator()(state_type state, ...)
  };
  or:
  class{
      void operator()(state_type && state, ...)
  };
 */
template <class T, class = void>
struct stepper_accepting_lvalue_state : std::true_type{};

template <class T>
struct stepper_accepting_lvalue_state<
  T,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >()
	)
       )
      >::value
    >
  > : std::false_type{};

/*
  for any stepper we need to ensure its operator() accepts the
  state by **non-const** reference. One (maybe only) way to detect
  that is check that if operator() can bind an rvalue state.
  Basically if we tried to bind an rvalue to an lvalue reference,
  it should fail, and if it fails that is good because it is what we want.

  In brief, we want the stepper to have this:
  class{
      void operator()(state_type & state, ...)
  };

  NOT this:
  class{
      void operator()(state_type state, ...)
  };
  or:
  class{
      void operator()(state_type && state, ...)
  };
 */

template <class T, class AuxT, class ...Args>
struct variadic_stepper_accepting_lvalue_state : std::true_type{};

template <class T, class AuxT, class ...Args>
struct variadic_stepper_accepting_lvalue_state<
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval<AuxT>(), std::declval<Args>()...
	)
       )
      >::value
    >, T, AuxT, Args...
  > : std::false_type{};



} // end namespace pressio::ode::impl


template <class T, class IndVarType, class StateType, class enable = void>
struct StateObserver : std::false_type{};

template <class T, class IndVarType, class StateType>
struct StateObserver<
  T, IndVarType, StateType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().operator()
	       (
		std::declval< ::pressio::ode::StepCount >(),
		std::declval< IndVarType >(),
		std::declval<StateType const &>()
		)
	       )
      >::value
    >
  > : std::true_type{};


template <class T, class IndVarType, class StateType, class enable = void>
struct StateGuesser : std::false_type{};

template <class T, class IndVarType, class StateType>
struct StateGuesser<
  T, IndVarType, StateType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>()
	       (
		std::declval< ::pressio::ode::StepCount >(),
		std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
		std::declval<StateType &>()
		)
	       )
      >::value
    && impl::guesser_taking_state_by_ref<T, IndVarType, StateType>::value
    >
  > : std::true_type{};


template <class T, class IndVarType, class Enable = void>
struct StepSizePolicy
  : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicy<
  T, IndVarType,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >()
	)
       )
      >::value
    && impl::step_size_policy_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};


template <class T, class IndVarType, class Enable = void>
struct StepSizePolicyWithReductionScheme
  : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicyWithReductionScheme<
  T, IndVarType,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeMinAllowedValue<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeScalingFactor<IndVarType> & >()
	)
       )
      >::value
    && impl::step_size_policy_with_reduc_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};

template <class T, class = void>
struct Steppable : std::false_type{};

template <class T>
struct Steppable<
  T,
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_independent_variable_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type & >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >()
	)
       )
      >::value
    && impl::stepper_accepting_lvalue_state<T>::value
    >
  > : std::true_type{};

template <class T> using StronglySteppable = Steppable<T>;


template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs : std::false_type{};

template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs<
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_independent_variable_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval<typename T::state_type & >(),
	std::declval<::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval<::pressio::ode::StepCount >(),
	std::declval<::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval< AuxT >(), std::declval<Args>()...
	)
       )
      >::value
    && impl::variadic_stepper_accepting_lvalue_state<void, T, AuxT, Args...>::value
    >,
  T, AuxT, Args...
  > : std::true_type{};

template <class T, class AuxT, class ...Args>
using StronglySteppableWithAuxiliaryArgs = SteppableWithAuxiliaryArgs<T, AuxT, Args...>;

}} // end namespace pressio::ode
#endif
