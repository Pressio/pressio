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
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (2)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser);

  template<
    class StepperType, class StateType, class GuesserType,
    class ObserverType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (3)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class ObserverType, class IndVarType
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (4)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType, class StateType, class GuesserType,
    class IndVarType, class AuxT, class ...Args
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (5)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType & guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class GuesserType, class IndVarType, class AuxT, class ...Args>
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (6)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType, class StateType, class GuesserType,
    class ObserverType, class IndVarType, class AuxT, class ...Args
    >
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (7)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
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
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer,
					     AuxT && auxArg,
					     Args && ... args);

  }} // end namespace pressio::ode

Parameters
----------

* ``stepper``: object that knows *how to* perform a single step.

* ``state``: object that represents the "state" to update

* ``startVal``: the independent variable starting value

* ``numSteps``: how many steps to take

* ``stepSizePolicy``: functor to set the step size

* ``stepSize``: *constant* step size to use for each step

* ``guesser``: functor to overwrite the state with a guess **before** doing a step

* ``observer``: object to "observe" the state's evolution, which can be used
  to potentially collect necessary data/metrics/statistics or do other things from the state.

* ``auxArg``, ``args``: arbitrary objects that are perfectly forwarded to the stepper's operator().


Constraints
-----------

* ``StepperType``:

  - for 1,2,3,5: must satisfy the `Steppable concept <ode_concepts/c6.html>`_

  - for 5,6,7,8, must satisfy the `SteppableWithAuxiliaryArgs concept <ode_concepts/c7.html>`_

* ``StepSizePolicyType`` must model the `StepSizePolicy concept <ode_concepts/c8.html>`_

* ``GuesserType`` must model the `StateGuesser concept <ode_concepts/c11.html>`_

* ``ObserverType`` must model he `StateObserver concept <ode_concepts/c10.html>`_

Preconditions
-------------

:red:`finish`

Mandates
--------

* ``std::is_same<IndVarType, typename StepperType::independent_variable_type>``

* ``std::is_same<StateType, typename StepperType::state_type>``

Return value
------------

None

Postconditions and Side Effects
-------------------------------

:red:`finish`
