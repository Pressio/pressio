.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_n_steps_with_pre_step_guesser``
=========================================

Header: ``<pressio/ode_advancers.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ode{

  template<
    class StepperType,
    class StateType,
    class GuesserType,
    class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
          && StateGuesser<GuesserType, IndVarType, StateType>
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (1)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser);

  template<
    class StepperType,
    class StateType,
    class StepSizePolicyType,
    class GuesserType,
    class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
	  && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
          && StateGuesser<GuesserType, IndVarType, StateType>
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (2)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser);

  template<
    class StepperType,
    class StateType,
    class GuesserType,
    class ObserverType,
    class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
          && StateGuesser<GuesserType, IndVarType, StateType>
	  && StateObserver<ObserverType &&, IndVarType, StateType>::value
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (3)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType,
    class StateType,
    class StepSizePolicyType,
    class GuesserType,
    class ObserverType,
    class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
	  && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
          && StateGuesser<GuesserType, IndVarType, StateType>
	  && StateObserver<ObserverType &&, IndVarType, StateType>::value
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (4)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer);

  template<
    class StepperType,
    class StateType,
    class GuesserType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
          && StateGuesser<GuesserType, IndVarType, StateType>
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (5)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType,
    class StateType,
    class StepSizePolicyType,
    class GuesserType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateGuesser<GuesserType, IndVarType, StateType>
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (6)
					     StateType & state,
					     const IndVarType & startVal,
					     StepSizePolicyType && stepSizePolicy,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType,
    class StateType,
    class GuesserType,
    class ObserverType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
          && StateGuesser<GuesserType, IndVarType, StateType>
          && StateObserver<ObserverType, IndVarType, StateType>
  #endif
  void advance_n_steps_with_pre_step_guesser(StepperType & stepper,         (7)
					     StateType & state,
					     const IndVarType & startVal,
					     const IndVarType & stepSize,
					     ::pressio::ode::StepCount numSteps,
					     GuesserType && guesser,
					     ObserverType && observer,
					     AuxT && auxArg,
					     Args && ... args);

  template<
    class StepperType,
    class StateType,
    class StepSizePolicyType,
    class GuesserType,
    class ObserverType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateGuesser<GuesserType, IndVarType, StateType>
	  && StateObserver<ObserverType, IndVarType, StateType>
  #endif
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

Description
-----------

Overload set to use a stepper object to update/advance a state :math:`n`
times such that *before* each step you can provide a guess for the solution.


Parameters
----------

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``stepper``
     - object that knows *how to* perform a single step

   * - ``state``
     - the "state" information to update

   * - ``startVal``
     - starting value of the independent variable

   * - ``numSteps``
     - how many steps to take

   * - ``stepSizePolicy``
     - functor to set the step size

   * - ``stepSize``
     - *constant* step size to use for each step

   * - ``guesser``
     - functor responsible to overwrite the state with a guess **before** doing a step

   * - ``observer``
     - object to "observe" the state's evolution, which can be used
       to potentially collect necessary data/metrics/statistics or do other things from the state

   * - ``auxArg``, ``args``
     - extra arguments that might be needed by the stepper to execute one step

Constraints
-----------

Concepts are documented `here <ode_concepts.html>`__.
Note: constraints are enforced via proper C++20 concepts when ``PRESSIO_ENABLE_CXX20`` is enabled,
otherwise via SFINAE and static asserts.

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
