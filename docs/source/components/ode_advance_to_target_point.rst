.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_to_target_point``
===========================

Header: ``<pressio/ode_advancers.hpp>``

API
---

.. code-block:: cpp

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
  #endif
   void advance_to_target_point(StepperType & stepper,          (1)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class ObserverType,
     class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires Steppable<StepperType>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
  #endif
   void advance_to_target_point(StepperType & stepper,          (2)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class IndVarType,
     class AuxT,
     class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>::value
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
	  && (!StateObserver<AuxT, IndVarType, StateType>)
  #endif
   void advance_to_target_point(StepperType & stepper,          (3)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       AuxT && auxArg,
			       Args && ... args);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class ObserverType,
     class IndVarType,
     class AuxT,
     class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires SteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>::value
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
	  && (!StateObserver<AuxT, IndVarType, StateType>)
  #endif

   void advance_to_target_point(StepperType & stepper,          (4)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer,
			       AuxT && auxArg,
			       Args && ... args);

Description
-----------

Overload set for using a stepper object to update/advance a state
until the independent variable reaches a target value.

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

   * - ``startVal``, ``finalVal``
     - the starting and final value of the independent variable

   * - ``numSteps``
     - how many steps to take

   * - ``stepSizePolicy``
     - functor to set the step size

   * - ``stepSize``
     - *constant* step size to use for each step

   * - ``observer``
     - object to "observe" the state's evolution, which can be used
       to potentially collect necessary data/metrics/statistics or do other things from the state

   * - ``auxArg``, ``args``
     - extra arguments that might be needed by the stepper to execute one step


Constraints
-----------

Each overload is associated with a set of constraints.
With C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:


* `Steppable <ode_concepts_various/steppable.html>`_

* `SteppableWithAuxiliaryArgs <ode_concepts_various/steppable_args.html>`_

* `StepSizePolicy <ode_concepts_various/step_size_pol.html>`_

* `StateObserver <ode_concepts_various/state_observer.html>`_

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
