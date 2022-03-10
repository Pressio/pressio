.. role:: raw-html-m2r(raw)
   :format: html

advancers
=========

.. note::

    Defined in header: ``<pressio/ode_advancers.hpp>``

    Public namespace: ``pressio::ode``

Scope and overview
------------------

.. note::

    Provides functionalities implementing various strategies
    to "update" some information (i.e., a state) via an object
    that satisfies a "steppable" concept (more details below).

Why is this useful? Suppose that you have a generic application/usecase
involving the concept of a discrete time or just discrete steps,
a notion of "state" that defines your information at a specific step,
and a "stepper" object, i.e. an object/abstraction that accepts a state
and knows how to take a "step" to update that state.
In such a scenario, it can be useful to have different strategies
to advance the state, for example: advancing for a fixed number of times,
advancing until a certain condition is met, etc.
This is what the ``pressio/ode_advancers`` provide.
We are growing this capability, so if you need something but it is
not alrady supported, open an issue or a PR on `github <https://github.com/Pressio/pressio>`_.

Obviously, this is not necessarily specific to applied math problems or
scientific computing applications, but can be something generic.
In the context of scientific computing, one might immediately recognize
these functions to be useful for *time integration* of dynamical systems.
This is the reason why these functionalities are currently inside ``pressio::ode``\ ,
and this is also why we specifically use "time" below.
Later on, we might move or generalize them.

API
---

.. note::

    Overload set for advancing for fixed number of steps

.. code-block:: cpp

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class ...Args
     >                                                                     (1)
   void advance_n_steps(StepperType & stepper,
                        StateType & state,
                        const TimeType start_time,
                        const TimeType time_step_size,
                        const pressio::ode::step_count_type num_steps,
                        Args && ... args);

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class StepSizeSetterType,
     class ...Args
     >                                                                     (2)
   void advance_n_steps(StepperType & stepper,
                        StateType & state,
                        const TimeType start_time,
                        StepSizeSetterType && time_step_size_manager,
                        const pressio::ode::step_count_type num_steps,
                        Args && ... args);

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class ObserverType,
     class ...Args
     >                                                                     (3)
   void advance_n_steps_and_observe(StepperType & stepper,
                                    StateType & state,
                                    const TimeType start_time,
                                    const TimeType time_step_size,
                                    const pressio::ode::step_count_type num_steps,
                                    ObserverType & observer,
                                    Args && ... args);

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class StepSizeSetterType,
     class ObserverType,
     class ...Args
     >                                                                     (4)
   void advance_n_steps_and_observe(StepperType & stepper,
                                    StateType & state,
                                    const TimeType start_time,
                                    StepSizeSetterType && time_step_size_manager,
                                    const pressio::ode::step_count_type num_steps,
                                    ObserverType & observer,
                                    Args && ... args);


* (1,2): overload set for advancing for a fixed number of steps
* (3,4): overload set for advancing for a fixed number of steps accepting
  also an "observer" to monitor the evolution of the state at each step (more on this below)

.. note::

    Overload set for advancing to target time

.. code-block:: cpp

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class StepSizeSetterType,
     class ...Args
     >
   void advance_to_target_time(StepperType & stepper,
                               StateType & state,
                               const TimeType start_time,
                               const TimeType final_time,
                               StepSizeSetterType && time_step_size_manager,
                               Args && ... args);

   template<
     class StepperType,
     class StateType,
     class TimeType,
     class StepSizeSetterType,
     class ObserverType,
     class ...Args
     >
   void advance_to_target_time_and_observe(StepperType & stepper,
                                           StateType & state,
                                           const TimeType start_time,
                                           const TimeType final_time,
                                           StepSizeSetterType && time_step_size_manager,
                                           ObserverType & observer,
                                           Args && ...args);

Parameters and Requirements
---------------------------

* 
  ``stepper``\ :

  * the steppable object, must conform to:

    .. code-block:: cpp

       class SteppableClass
       {
       public:
       void operator()(StateType & state,
                       const TimeType current_time,
                       const TimeType time_step_size_to_use,
                       const int32_t step_count
                       /*[, other arguments needed: these will be forwarded to your
                       stepper object by the advance methods as shown above ]*/)
       {
        // do what you need to update state
       }
       };

* 
  ``state``\ : self-explanatory

* 
  ``start_time``\ : self-explanatory

* 
  ``final_time``\ : self-explanatory

* 
  ``num_steps``\ : self-explanatory

* 
  ``time_ste_size_manager``\ :

  * functor for setting the time step size

    .. code-block:: cpp

       class StepSizeSetter
       {
       public:
       void operator()(const int32_t step_count,
                       const TimeType time,
                       TimeType & dt) const
       {
        // update somehow dt
       }
       };

* 
  ``observer``\ :

  * functor to "observe" the state during the time integration,
    allowing you to potentially collect necessary data/metrics/statistics.
    Must conform to:

    .. code-block:: cpp

       class ObserverClass
       {
       public:
       void operator()(const int32_t & step_count,
                       const TimeType current_time,
                       const StateType & current_state) const;
       };

* 
  ``args...``\ : optional objects to forward to stepper's ``operator()`` that are
  potentially needed to perform one step.
  Note that these are optional, because your stepper might not need anything.
  The advance functions will simply forward all these to the ``operator()``
  of the stepper object. Is having these optional argument really needed?
  One might argue that if the stepper needs to have access to
  some objects, then these auxiliary objects can be passed to the stepper
  upon construction, and so the stepper would already have access to them
  when its ``operator()`` is called.
  This is true *if* there is a valid reason for making these auxiliary object
  data members of the steppable class.
  The problem with this approach is that the stepper would need to know
  about all theses other types, and this creates a very tight
  coupling between different objects,
  that might be not even needed.
  By letting users provide a pack to the advance functions,
  we are able to decouple this structure such that any object
  *only* needed within the stepper's ``operator()`` can be directly passed.
  To explain this, look at the example below.

When are the variadic arguments useful?
---------------------------------------

For example in a scenario like the following:

.. code-block:: cpp

   struct SteppableClass
   {
     template<class AuxiliaryType>
     void operator()(StateType & state,
                     const TimeType current_time,
                     const TimeType time_step_size_to_use,
                     const int32_t step_count,
                     AuxiliaryType & aux)
     {
        const int value = aux.doComplicateCalculation();
        if (value % 2 == 0){
          // update state somehow
        }
        else{
          // update state in a different way
        }
     }
   };

   class Foo{
     int doComplicateCalculation(){}
   };

   class Bar{
     int doComplicateCalculation(){}
   };

   template<class AuxType>
   void run()
   {
     AuxType a;
     advance_n_steps(stepper, /*state*/, /*start t*/, /*step size*/, /*num steps*/, a);
   }

   int main()
   {
     SteppableClass stepper;
     run<Foo>();
     run<Bar>();
   }

Here, we don't want to parametrize the ``StepperClass``
on the auxiliary class needed because we don't need it as member.
We just need to access the auxiliary class inside the ``operator()``.
