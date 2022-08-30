.. role:: raw-html-m2r(raw)
   :format: html

``advance_to_target_point_with_step_recovery``
==============================================

Defined in header: ``<pressio/ode_advancers.hpp>``

Scope
-----

Overload set for using a stepper object to update/advance a state
until the independent variable reaches a target value with support
for recovering if the step fails.

API
---

.. code-block:: cpp

   template<class StepperType, class StateType, class StepSizePolicyType, class IndVarType>
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (1)
						  StateType & state,
						  const IndVarType & start_val,
						  const IndVarType & final_val,
						  StepSizePolicyType && step_size_policy);

   template<
     class StepperType, class StateType, class StepSizePolicyType,
     class ObserverType, class IndVarType>
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (2)
						  StateType & state,
						  const IndVarType & start_val,
						  const IndVarType & final_val,
						  StepSizePolicyType && step_size_policy,
						  ObserverType && observer);

   template<
     class StepperType, class StateType, class StepSizePolicyType,
     class IndVarType, class AuxT, class ...Args>
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (3)
						  StateType & state,
						  const IndVarType & start_val,
						  const IndVarType & final_val,
						  StepSizePolicyType && step_size_policy,
						  AuxT && auxArg,
						  Args && ... args);

   template<
     class StepperType, class StateType, class StepSizePolicyType,
     class ObserverType, class IndVarType, class AuxT, class ...Args>
   void advance_to_target_point_with_step_recovery(StepperType & stepper,          (4)
						  StateType & state,
						  const IndVarType & start_val,
						  const IndVarType & final_val,
						  StepSizePolicyType && step_size_policy,
						  ObserverType && observer,
						  AuxT && auxArg,
						  Args && ... args);


Parameters and Requirements
---------------------------

* ``stepper``: an object that knows *how to* perform a single step.

  - for overloads 1,2 ``StrongStepperType`` must satisfy the `StronglySteppable concept <ode_concepts/c6.html#stronglysteppable>`_

  - for overloads 3,4 ``StrongStepperType`` must satisfy the `StronglySteppableWithAuxiliaryArgs concept <ode_concepts/c7.html#stronglysteppablewithauxiliaryargs>`_

* ``state``: self-explanatory

  - constraint: ``std::is_same<StateType, typename StepperType::state_type>``

* ``start_val``, ``final_val``: self-explanatory

  - constraint: ``std::is_same<IndVarType, typename StepperType::independent_variable_type>``

* ``step_size_policy``: functor to set the step size

  - must conform to the `StepSizePolicyWithReductionScheme concept <ode_concepts/c9.html>`_

* ``observer``: functor to "observe" the state's evolution

  - must conform to the `StateObserver concept <ode_concepts/c10.html>`_

  - scope: to potentially collect necessary data/metrics/statistics or
    do other things from the state.

* ``auxArg``, ``args``: arbitrary objects that are perfectly forwarded to the stepper's operator().
