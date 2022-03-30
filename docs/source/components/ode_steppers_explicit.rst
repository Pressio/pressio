.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

explicit steppers
=================

.. admonition:: Info
   :class: important

   Header: ``<pressio/ode_steppers_explicit.hpp>``

   Public namespace: ``pressio::ode``


.. admonition:: Description

   :medium:`Classes representing different methods
   to execute a *single step* when doing explicit time integration.`

Overview
--------

Suppose that you have a system of the form:

.. math::

    \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...)

where :math:`y` is the state, :math:`f` is the RHS (also called velocity below), :math:`t` is time.\ :raw-html-m2r:`<br/>`
Explicit methods calculate the state of a system at a later time
from the state at the current time and potentially previous times.
A pressio "explicit stepper" is an abstraction that represents the "how" to take a step.

API, Parameters and Requirements
--------------------------------

.. code-block:: cpp

   namespace pressio{ namespace ode{

   template<class StateType, class SystemType>
   auto create_explicit_stepper(pressio::ode::StepScheme scheme_name,
                                const StateType & state,
                                const SystemType & system);

   }}

*
  ``scheme_name``

  * one of the following enum values ``pressio::ode::StepScheme::{ForwardEuler, RungeKutta4, AdamsBashforth2, SSPRungeKutta3}``.

*
  ``state``\ :

  * an instance of your state. The type must be copy constructible.

*
  ``system``\ :

  *
    an instance of your problem class that encapsulated how to compute the velocity :math:`f` for a given state and time.
    The system class must conform to the following API:

    .. code-block:: cpp

        struct SystemForExplicitStepper
        {
	  using scalar_type   = /* */;
	  using state_type    = /* */;
	  using velocity_type = /* */;

	  velocity_type createVelocity() const;
	  void velocity(const state_type &, scalar_type time, velocity_type &) const;
        };

  the nested aliases ``scalar_type``\ , ``state_type`` and ``velocity_type`` must be *valid* types since
  they are detected by pressio

*
  if ``StateType`` is the type deduced for ``state`` from ``create_...``\ , the following must hold:\ :raw-html-m2r:`<br/>`
  ``std::is_same<StateType, typename SystemForExplicitOde::state_type>::value == true``



Stepper Class Public API
------------------------

A stepper class meets the following API:

.. code-block:: cpp

    template<class ScalarType, class StateType, class ... Args>
    class ExplicitStepper
    {
       template<class StepCountType>
       void operator()(StateType odeState,
		       const ScalarType & currentTime,
		       const ScalarType & dt,
		       StepCountType stepNumber);
    }


.. tip::

   The stepper class satisfies the "steppable" concept discussed `here <ode_advance.html>`_\ , so one can use the "advancers" functions to step forward.


Example code
------------

.. code-block:: cpp

   #include "pressio/type_traits.hpp"
   #include "pressio/ode_advancers.hpp"
   #include "pressio/ode_steppers_explicit.hpp"
   int main()
   {
     // assuming that:
     // myState  is the state
     // mySystem is the system instance

     namespace pode = pressio::ode;
     const auto scheme = pode::StepScheme::ForwardEuler;
     auto stepper = pode::create_explicit_stepper(scheme, myState, mySystem);

     // use the stepper to advance in time,
     // for example using the advancer function
     const double time0 = 0.;
     const double dt = 0.1;
     const pode::step_count_type num_steps = 100;
     pode::advance_n_steps(stepper, myState, time0, dt, num_steps);
   }


|


Required specializations for custom types
-----------------------------------------

When using custom data types not supported in `pressio ops <ops.html>`_\ , you need to provide specializations of a trait class and certain operations
and make them "visible" to the compiler to find them and such that pressio can operate on your data.
For the sake of explanation, suppose that you use ``double``
as value type and ``ACustomStateType`` is what you use for the state, then you would need to do something like this:

.. code-block:: cpp

   #include "pressio/type_traits.hpp"

   // assuming ACustomStateType has already been declared

   namespace pressio{

   template<> struct Traits<ACustomStateType>{
     using scalar_type = double;
   };

   namespace ops{

   void deep_copy(ACustomStateType & dest, const ACustomStateType & src){
     /* deep copy src into dest */
   }

   ACustomStateType clone(const ACustomStateType & src){
     /* return a deep copy of src */
   }

   void set_zero(ACustomStateType & object){
     /* set elements to zero */
   }

   void update(ACustomStateType & v,        const double a,
               const ACustomStateType & v1, const double b)
   {
     // elementwise compute : v = a*v + b*v1
   }

   void update(ACustomStateType & v,        const double a,
               const ACustomStateType & v1, const double b,
               const ACustomStateType & v2, const double c)
   {
     // elementwise compute : v = a*v + b*v1 + c*v2
   }

   void update(ACustomStateType & v,        const double a,
               const ACustomStateType & v1, const double b,
               const ACustomStateType & v2, const double c,
               const ACustomStateType & v3, const double d)
   {
     // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3
   }

   void update(ACustomStateType & v,        const double a,
               const ACustomStateType & v1, const double b,
               const ACustomStateType & v2, const double c,
               const ACustomStateType & v3, const double d,
               const ACustomStateType & v4, const double e)
   {
     // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3 + e*v4
   }
   }}//end namepsace pressio::ops

   #include "pressio/ode_advancers.hpp"
   #include "pressio/ode_steppers_explicit.hpp"

   int main()
   {
     // same code as shown above
   }

Note that in the snippet above the order of the include statements matter:
this is because your ``Trait`` and kernel specializations need to be found by the compiler.
However, to make the code cleaner, you can obviously move all kernels specializations
to a separate header file, but make sure to keep the correct order, for example as follows:

.. code-block:: cpp

   #include "pressio/type_traits.hpp"
   #include "my_specializations.hpp" // contains all your specializations
   #include "pressio/ode_advancers.hpp"
   #include "pressio/ode_steppers_explicit.hpp"
   int main()
   {
     // same code as shown above
   }
