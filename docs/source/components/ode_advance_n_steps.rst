.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_n_steps``
===================

Header: ``<pressio/ode_advancers.hpp>``

Scope
-----

Overload set for using a stepper object to update/advance a state :math:`n` times.

API
---

.. code-block:: cpp

  namespace pressio { namespace ode{

  template<class StepperType, class StateType, class IndVarType>
  void advance_n_steps(StepperType & stepper,                        (1)
                       StateType & state,
                       const IndVarType & startVal,
                       const IndVarType & stepSize,
                       const ::pressio::ode::StepCount numSteps);

  template<class StepperType, class StateType, class IndVarType, class StepSizePolicy>
  void advance_n_steps(StepperType & stepper,                        (2)
                       StateType & state,
                       const IndVarType & startVal,
                       const StepSizePolicy & stepSizePolicy,
                       const ::pressio::ode::StepCount numSteps);

  template<class StepperType, class StateType, class IndVarType, class ObserverType>
  void advance_n_steps(StepperType & stepper,                        (3)
                       StateType & state,
                       const IndVarType & startVal,
                       const IndVarType & stepSize,
                       const ::pressio::ode::StepCount numSteps,
		       ObserverType && observer);

  template<
    class StepperType, class StateType, class IndVarType,
    class StepSizePolicy, class ObserverType>
  void advance_n_steps(StepperType & stepper,                        (4)
                       StateType & state,
                       const IndVarType & startVal,
                       const StepSizePolicy & stepSizePolicy,
                       const ::pressio::ode::StepCount numSteps,
		       ObserverType && observer);

  template<
    class StepperType, class StateType, class IndVarType,
    class AuxT, class ...Args>
  void advance_n_steps(StepperType & stepper,                        (5)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps,
		       AuxT && auxArg,
		       Args && ... args);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class IndVarType, class AuxT, class ...Args>
  void advance_n_steps(StepperType & stepper,                        (6)
		       StateType & state,
		       const IndVarType & startVal,
		       StepSizePolicyType && stepSizePolicy,
		       ::pressio::ode::StepCount numSteps,
		       AuxT && auxArg,
		       Args && ... args);

  template<
    class StepperType, class StateType, class ObserverType,
    class IndVarType, class AuxT, class ...Args>
  void advance_n_steps(StepperType & stepper,                        (7)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps,
		       ObserverType && observer,
		       AuxT && auxArg,
		       Args && ... args);

  template<
    class StepperType, class StateType, class StepSizePolicyType,
    class ObserverType, class IndVarType, class AuxT, class ...Args>
  void advance_n_steps(StepperType & stepper,                        (8)
		       StateType & state,
		       const IndVarType & startVal,
		       StepSizePolicyType && stepSizePolicy,
		       ::pressio::ode::StepCount numSteps,
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

* ``observer``: object to "observe" the state's evolution, which can be used
  to potentially collect necessary data/metrics/statistics or do other things from the state.

* ``auxArg``, ``args``: arbitrary objects that are perfectly forwarded to the stepper's operator().


Constraints
-----------

* ``StepperType``:

  - for 1,2,3,5: must satisfy the `Steppable concept <ode_concepts/c6.html>`_

  - for 5,6,7,8, must satisfy the `SteppableWithAuxiliaryArgs concept <ode_concepts/c7.html>`_

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









.. *
..   ``args...``\ : optional objects to forward to stepper's ``operator()`` that are
..   potentially needed to perform one step.
..   Note that these are optional, because your stepper might not need anything.
..   The advance functions will simply forward all these to the ``operator()``
..   of the stepper object. Is having these optional argument really needed?
..   One might argue that if the stepper needs to have access to
..   some objects, then these auxiliary objects can be passed to the stepper
..   upon construction, and so the stepper would already have access to them
..   when its ``operator()`` is called.
..   This is true *if* there is a valid reason for making these auxiliary object
..   data members of the steppable class.
..   The problem with this approach is that the stepper would need to know
..   about all theses other types, and this creates a very tight
..   coupling between different objects,
..   that might be not even needed.
..   By letting users provide a pack to the advance functions,
..   we are able to decouple this structure such that any object
..   *only* needed within the stepper's ``operator()`` can be directly passed.
..   To explain this, look at the example below.

.. |

.. When are the variadic arguments useful?
.. ---------------------------------------

.. For example in a scenario like the following:

.. .. code-block:: cpp

..    struct SteppableClass
..    {
..      template<class AuxiliaryType>
..      void operator()(StateType & state,
..                      const TimeType current_time,
..                      const TimeType time_step_size_to_use,
..                      const int32_t step_count,
..                      AuxiliaryType & aux)
..      {
..         const int value = aux.doComplicateCalculation();
..         if (value % 2 == 0){
..           // update state somehow
..         }
..         else{
..           // update state in a different way
..         }
..      }
..    };

..    class Foo{
..      int doComplicateCalculation(){}
..    };

..    class Bar{
..      int doComplicateCalculation(){}
..    };

..    template<class AuxType>
..    void run()
..    {
..      AuxType a;
..      advance_n_steps(stepper, /*state*/, /*start t*/, /*step size*/, /*num steps*/, a);
..    }

..    int main()
..    {
..      SteppableClass stepper;
..      run<Foo>();
..      run<Bar>();
..    }

.. Here, we don't want to parametrize the ``StepperClass``
.. on the auxiliary class needed because we don't need it as member.
.. We just need to access the auxiliary class inside the ``operator()``.
