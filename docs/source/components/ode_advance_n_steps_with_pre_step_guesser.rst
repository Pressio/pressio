.. role:: raw-html-m2r(raw)
   :format: html

``advance_n_steps_with_pre_step_guesser``
=========================================

Header: ``<pressio/ode_advancers.hpp>``

Scope
-----

Overload set to use a stepper object to update/advance a state :math:`n`
times such that *before* each step you can provide a guess for the solution.

API
---

.. code-block:: cpp

  namespace pressio { namespace ode{

  template<
    class StepperType, class StateType, class GuesserType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (1)
					     StateType & state,
					     const IndVarType & start_val,
					     const IndVarType & step_size,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (2)
					     StateType & state,
					     const IndVarType & start_val,
					     StepSizePolicyType && step_size_policy,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser);

  template<
    class StepperType, class StateType, class GuesserType,
    class ObserverType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (3)
					     StateType & state,
					     const IndVarType & start_val,
					     const IndVarType & step_size,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class ObserverType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (4)
					     StateType & state,
					     const IndVarType & start_val,
					     StepSizePolicyType && step_size_policy,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType, class StateType, class GuesserType,
    class IndVarType, class AuxT, class ...Args
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (5)
					     StateType & state,
					     const IndVarType & start_val,
					     const IndVarType & step_size,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType & guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class IndVarType, class AuxT, class ...Args>
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (6)
					     StateType & state,
					     const IndVarType & start_val,
					     StepSizePolicyType && step_size_policy,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType, class StateType, class GuesserType,
    class ObserverType, class IndVarType, class AuxT, class ...Args
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (7)
					     StateType & state,
					     const IndVarType & start_val,
					     const IndVarType & step_size,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType & guesser,
					     ObserverType && observer,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class ObserverType, class IndVarType,
    class AuxT, class ...Args>
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (8)
					     StateType & state,
					     const IndVarType & start_val,
					     StepSizePolicyType && step_size_policy,
					     ::pressio::ode::StepCount num_steps,
					     GuesserType && guesser,
					     ObserverType && observer,
					     AuxT && auxArg,
					     Args && ... args);

  }} // end namespace pressio::ode


Parameters and Requirements
---------------------------

* ``stepper``: an object that knows *how to* perform a single step.

  - for overloads 1,2,3,4, ``StepperType`` must satisfy the `Steppable concept <ode_concepts/c6.html>`_

  - for overloads 5,6,7,8, ``StepperType`` must satisfy the `SteppableWithAuxiliaryArgs concept <ode_concepts/c7.html>`_

* ``state``: self-explanatory

  - constraint: ``std::is_same<StateType, typename StepperType::state_type>``

* ``start_val``: self-explanatory

  - constraint: ``std::is_same<IndVarType, typename StepperType::independent_variable_type>``

* ``num_steps``: self-explanatory

* ``step_size_policy``: functor to set the step size

  - must conform to the `StepSizePolicy concept <ode_concepts/c8.html>`_

* ``step_size``: the *constant* step size to use for each step

* ``observer``: functor to "observe" the state's evolution

  - must conform to the `StateObserver concept <ode_concepts/c10.html>`_

  - scope: to potentially collect necessary data/metrics/statistics or
    do other things from the state.

* ``guesser``: functor to overwrite the state with a guess **before** doing a step.

  - must conform to the `StateGuesser concept <ode_concepts/c11.html>`_

* ``auxArg``, ``args``: arbitrary objects that are perfectly forwarded to the stepper's operator().
