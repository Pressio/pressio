.. role:: raw-html-m2r(raw)
   :format: html

``advance_to_target_point``
===========================

Header: ``<pressio/ode_advancers.hpp>``


Scope
-----

Overload set for using a stepper object to update/advance a state
until the independent variable reaches a target value.

API
---

.. code-block:: cpp

   template<
     class StepperType, class StateType, class StepSizePolicyType, class IndVarType
     >
   void advance_to_target_point(StepperType & stepper,          (1)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & final_val,
			       StepSizePolicyType && stepSizePolicy);

   template<
     class StepperType, class StateType, class StepSizePolicyType,
     class ObserverType, class IndVarType>
   void advance_to_target_point(StepperType & stepper,          (2)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & final_val,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer);

   template<
     class StepperType, class StateType, class StepSizePolicyType, class IndVarType,
     class AuxT, class ...Args>
   void advance_to_target_point(StepperType & stepper,          (3)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & final_val,
			       StepSizePolicyType && stepSizePolicy,
			       AuxT && auxArg,
			       Args && ... args);

   template<
     class StepperType, class StateType, class StepSizePolicyType,
     class ObserverType, class IndVarType,
     class AuxT, class ...Args>
   void advance_to_target_point(StepperType & stepper,          (4)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & final_val,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer,
			       AuxT && auxArg,
			       Args && ... args);

Parameters
----------

* ``stepper``: object that knows *how to* perform a single step.

* ``state``: object that represents the "state" to update

* ``startVal``, ``final_val``: the independent variable starting value and target final value

* ``stepSizePolicy``: functor to set the step size

* ``observer``: object to "observe" the state's evolution, which can be used
  to potentially collect necessary data/metrics/statistics or do other things from the state.

* ``auxArg``, ``args``: arbitrary objects that are perfectly forwarded to the stepper's operator().

Constraints
-----------

* ``StepperType``:

  - for 1,2: must satisfy the `Steppable concept <ode_concepts/c6.html>`_

  - for 3,4 must satisfy the `SteppableWithAuxiliaryArgs concept <ode_concepts/c7.html>`_

* ``StepSizePolicyType`` must model the `StepSizePolicy concept <ode_concepts/c8.html>`_

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
