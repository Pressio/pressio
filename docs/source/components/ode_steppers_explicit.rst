.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

Explicit steppers
=================

Header: ``<pressio/ode_steppers_explicit.hpp>``

Public namespace: ``pressio::ode``

Scope
-----

A pressio "explicit stepper" is an abstraction representing
"how" to take a step when applying an :orange:`explicit scheme`
to initial value problems expressable in the following form

.. math::
  :label: ode_explicit_system1

    \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0

or with a mass matrix:

.. math::
  :label: ode_explicit_system2

    M(\boldsymbol{y}, t, ...) \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0


where :math:`y` is the state, :math:`f` is the right hand side (RHS),
:math:`t` is the independent variable,
and :math:`M` is the mass matrix. Note that both the right hand side
and the mass matrix potentially depend on the state and the independent variable.
We adopt the typical notation :math:`t` for the
independent variable but this does NOT necessarily mean that it refers to
time. The independent variable can represent something else.

.. admonition:: :medium:`Recall the definition`

   :medium:`Explicit methods calculate the next state using the current state and potentially previous ones.`


API
---

.. code-block:: cpp

   namespace pressio{ namespace ode{

   template<class SystemType>
   /*impl defined*/ create_explicit_stepper(StepScheme odeScheme,                   (1)
                                            SystemType && system);

   template<class SystemType, class MassMatrixOperatorType>
   /*impl defined*/ create_explicit_stepper(StepScheme odeScheme,                   (2)
                                            SystemType && system,
					    MassMatrixOperatorType && massMatOperator);

   }} //end namespace pressio::ode

Parameters
~~~~~~~~~~

* ``odeScheme``: the target integration scheme

  * choices: ``pressio::ode::StepScheme::{ForwardEuler, RungeKutta4, AdamsBashforth2, SSPRungeKutta3}``.

* ``system``: an object to compute the RHS for your problem

* ``massMatOperator``: the mass matrix operator

Constraints
~~~~~~~~~~~

- ``mpl::remove_cvref_t<SystemType>``: must model the `SystemWithRhs concept <ode_concepts/c1.html>`__
  or one subsuming it

- ``mpl::remove_cvref_t<MassMatrixOperatorType>``: must model either the `MassMatrixOperator concept <ode_concepts/c3a.html>`__
  or the `ConstantMassMatrixOperator concept <ode_concepts/c3b.html>`__

Preconditions
~~~~~~~~~~~~~

- if ``system`` does *not* bind to a temporary object,
  it must bind to an lvalue object whose lifetime is
  *longer* that that of the instantiated stepper, i.e., it is destructed
  *after* the stepper goes out of scope

- if passing an lvalue to ``massMatOperator``, the same precondition
  above applies

Mandates
~~~~~~~~

- if ``MassMatrixOperatorType`` meets the `MassMatrixOperator concept <ode_concepts/c3a.html>`__,
  then these must hold:

  - ``std::is_same< typename mpl::remove_cvref_t<SystemType>::independent_variable_type,
    typename mpl::remove_cvref_t<MassMatrixOperatorType>::independent_variable_type >::value == true``

  - ``std::is_same< typename mpl::remove_cvref_t<SystemType>::state_type,
    typename mpl::remove_cvref_t<MassMatrixOperatorType>::state_type >::value == true``

Return value
~~~~~~~~~~~~

An instance of a pressio explicit stepper suitable for the target scheme.
The return type is implementation defined.


Postconditions and side effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on what kind of problem (with or without the mass matrix)
you are trying to solve, the stepper object returned by pressio exposes different methods.

- If your system is of the form :eq:`ode_explicit_system1` and you use overload (1),
  the stepper instantiated by pressio exposes the following API:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class ExplicitStepper
      {
        public:
	  using state_type                = /* same type of your problem */;
	  using independent_variable_type = /* same as your problem */;

	  void operator()(StateType & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<IndVarType> & /**/);
      };


- If your system is of the form :eq:`ode_explicit_system2` and you use overload (2),
  and the stepper instantiated exposes the following API:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class ExplicitStepper
      {
        public:
	  using state_type                = /* same type of your problem */;
	  using independent_variable_type = /* same as your problem */;

	  template<class LinearSolverType>
	  void operator()(StateType & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<IndVarType> & /**/,
			  LinearSolverType & solver);
      };

  Note that, if you have a mass matrix, the stepper needs a linear solver
  to execute one step. The linear solver should conform to [this API](:red:`TODO`).
  Why do we need a linear solver? This is for performance reasons. One possibility would
  be to directly invert the mass matrix and reformulate the original problem as :math:`dy/dt = M^{-1} f`.
  This is impractical to do for most real applications where the mass matrix is large.
  A better approach, which avoids the inverse, is to use the mass matrix to solve
  sequences of problems of the form :math:`M x = f` and then use this :math:`x` in a certain way
  to update the state :math:`y`. Obviously, the details of how this is done depend on
  the chosen scheme. You are just responsible of passing an object to do a linear solve
  of such system. Since you know what matrix you have, its structure and what
  is the right hand side, you can then provide the most suitable linear solver.


- if the function is called with an rvalue for either of the arguments,
  the constructor of the stepper
  will try to use move semantics. If move semantics are implemented, the temporary
  is moved from and no new memory allocation occurs. If move semantics fallback to copying,
  then a copy of the original object is made. Either way, the stepper instantiated
  directly manages the lifetime of the object.

- if the function is called with an lvalue for either of the arguments,
  the constructor of the stepper
  will create a *reference to the* object. In such case, the stepper
  object does NOT manage the lifetime of it. You are responsible to
  ensure that the ``system`` or the ``massMatOperator`` are destructed
  *after* the stepper goes out of scope

Use the stepper
---------------

Now that we discussed what a stepper is and what
methods it exposes, we have everything we need to
discuss how to use it to advance in time.
For example, you can either implement your step
loop which sequentially calls the stepper.
However, the key thing to notice here is that a stepper
satisfies the "steppable" concept
discussed `here <ode_concepts/c6.html>`_\ , so one can
use the "advancers" functions to step forward.

Let's first look at an example for a system without mass matrix:

.. code-block:: cpp

   #include "pressio/type_traits.hpp"
   #include "pressio/ode_advancers.hpp"
   #include "pressio/ode_steppers_explicit.hpp"
   int main()
   {
     MyRhsEvalClass myRhsEval(/*whatever args*/);

     namespace pode = pressio::ode;
     const auto scheme = pode::StepScheme::ForwardEuler;
     auto stepper = pode::create_explicit_stepper(scheme, myRhsEval);

     auto myState = myRhsEval.createState();
     // initialize myState to initial condition

     using time_type = typename MyRhsEvalClass::independent_variable_type;
     const time_type t0 = 0.;
     const time_type dt = 0.1;
     const int num_steps = 100;
     pode::advance_n_steps(stepper, myState, t0, dt, num_steps);
   }


Let's look at an example for a system *with* mass matrix:

.. code-block:: cpp

   #include "pressio/type_traits.hpp"
   #include "pressio/ode_advancers.hpp"
   #include "pressio/ode_steppers_explicit.hpp"
   int main()
   {
     MyRhsEvalClass myRhsEval(/*whatever args*/);
     MyMassMatrixOpClass myMassMatrixOp(/*whatever args*/);

     namespace pode = pressio::ode;
     const auto scheme = pode::StepScheme::ForwardEuler;
     auto stepper = pode::create_explicit_stepper(scheme, myRhsEval, myMassMatrixOp);

     auto myState = myRhsEval.createState();
     // initialize myState to initial condition

     auto linearSolver = /* create linear solver */;

     using time_type = typename MyRhsEvalClass::independent_variable_type;
     const time_type time0 = 0.;
     const time_type dt = 0.1;
     const int num_steps = 100;
     pode::advance_n_steps(stepper, myState, time0, dt, num_steps, linearSolver);
   }
