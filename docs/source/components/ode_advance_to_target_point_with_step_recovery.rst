.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_to_target_point_with_step_recovery``
==============================================

Defined in header: ``<pressio/ode_advancers.hpp>``

API
---

.. code-block:: cpp

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class IndVarType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires StronglySteppable<StepperType>
	   && StepSizePolicyWithReductionScheme<StepSizePolicyType, IndVarType>
   #endif
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (1)
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
     requires StronglySteppable<StepperType>
	   && StepSizePolicyWithReductionScheme<StepSizePolicyType, IndVarType>
	   && StateObserver<ObserverType&&, IndVarType, StateType>::value
   #endif
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (2)
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
     requires StronglySteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
           && StepSizePolicyWithReductionScheme<StepSizePolicyType, IndVarType>
	   && (!StateObserver<AuxT, IndVarType, StateType>)
   #endif
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (3)
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
     requires StronglySteppableWithAuxiliaryArgs<StepperType, AuxT, Args...>
           && StepSizePolicyWithReductionScheme<StepSizePolicyType, IndVarType>
	   && StateObserver<ObserverType, IndVarType, StateType>
   void
   #endif
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (4)
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
until the independent variable reaches a target value with support
for recovering if the step fails.


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

With C++20, the constraints would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message) and/or SFINAE.

The concepts are documented `here <ode_concepts.html>`__.

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
