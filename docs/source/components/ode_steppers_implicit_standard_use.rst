.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

Standard Use
============

Header: ``<pressio/ode_steppers_implicit.hpp>``

Public namespace: ``pressio::ode``


Scope
-----

This page describes the "standard use" of an "implicit stepper"
for systems expressed either as

.. math::
  :label: ode_implicit_system_std1

    \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0

or with a mass matrix:

.. math::
   :label: ode_implicit_system_std2

    M(\boldsymbol{y}, t, ...) \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0


What do we mean by `standard` use?
----------------------------------

We refer to a use case as "standard" when you define your problem
via objects that are responsible of computing ``f``, its jacobian ``df/dy``
and, if applicable, the mass matrix ``M``,
and you want to do implicit stepping for it.

API
---

.. code-block:: cpp

   namespace pressio{ namespace ode{

   template<class SystemType>
   /*impl defined*/ create_implicit_stepper(StepScheme odeScheme,
                                            const SystemType & system);

   template<class SystemType, class MassMatrixOperatorType>
   /*impl defined*/ create_implicit_stepper(StepScheme odeScheme,
                                            const SystemType & system,
					    MassMatrixOperatorType && massMatOperator);

   }} //end namespace pressio::ode


Parameters
~~~~~~~~~~

* ``odeScheme``: the target integration scheme

  * choices: ``pressio::ode::StepScheme::{BDF1, BDF2, CrankNicolson}``

* ``system``: an object that evaluates the rhs and its jacobian for your problem

* ``massMatOperator``: the mass matrix operator


Constraints
~~~~~~~~~~~

- ``SystemType``: must model the `SystemWithRhsAndJacobian concept <ode_concepts/c2.html>`__
  or one subsuming it

- ``std::decay_t<MassMatrixOperatorType>``: must model either the `MassMatrixOperator concept <ode_concepts/c3a.html>`__
  or the `ConstantMassMatrixOperator concept <ode_concepts/c3b.html>`__


Preconditions
~~~~~~~~~~~~~

- ``system`` must bind to an lvalue object whose lifetime is
  longer that that of the instantiated stepper, i.e., it is destructed
  *after* the stepper goes out of scope

- if passing an lvalue reference for ``massMatOperator``, the same precondition
  above applies

Mandates
~~~~~~~~

- if ``MassMatrixOperatorType`` meets the `MassMatrixOperator concept <ode_concepts/c3a.html>`__,
  then the following must hold:

  - ``std::is_same< typename SystemType::independent_variable_type,
    typename std::decay_t<MassMatrixOperatorType>::independent_variable_type >::value == true``

  - ``std::is_same< typename SystemType::state_type,
    typename std::decay_t<MassMatrixOperatorType>::state_type >::value == true``

Return value
~~~~~~~~~~~~

An instance of a pressio implicit stepper suitable for the target scheme.
The return type is implementation defined.

Postconditions and side effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The returned stepper object is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class ImplicitStepper
    {
      public:
        // these nested type aliases must be here
        using independent_variable_type = /* same as in your system class */;
        using state_type                = /* same as in your system class */;
        using residual_type             = /* equal to right_hand_side_type in your system class */;
        using jacobian_type	        = /* same as in your system class */;

        template<class SolverType>
        void operator()(StateType & /**/,
			const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			pressio::ode::StepCount /**/,
			pressio::ode::StepSize<independent_variable_type> /**/,
			SolverType & /**/)

        residual_type createResidual() const;
        jacobian_type createJacobian() const;
        void residualAndJacobian(const state_type & odeState,
	                         residual_type & R,
				 jacobian_type & J,
				 bool computeJacobian) const;
    };


- if ``create_implicit_stepper()`` was called with an
  rvalue for the mass matrix operator, the constructor of the stepper object
  will try to use move semantics. If move semantics are implemented, the temporary
  is moved from and no new memory allocation occurs. If move semantics fallback to copying,
  then a copy of the original object is made. Either way, the stepper instantiated
  directly manages the lifetime of the object.

- if ``create_implicit_stepper()`` was called with an
  lvalue for the mass matrix operator, the stepper will create a *reference* to this object.
  In such case, the stepper instance does NOT manage the lifetime of it.
  You are responsible to ensure that the ``massMatOperator`` is destructed
  *after* the stepper goes out of scope

- the returned stepper instance holds a reference to ``system``,
  therefore the stepper does NOT manage its lifetime. You need to ensure
  that the ``system`` is destructed *after* the stepper goes out of scope


Use the stepper
---------------

The stepper created using the functions satisfies two concepts:

- the ``SteppableWithAuxiliaryArgs`` concept discussed `here <ode_concepts/c7.html>`__.

- the ``SystemWithFusedResidualAndJacobian`` concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.

Therefore, to advance the stepper you can use a solver from the nonlinear_solvers component and
the "advance" functions to step forward. Or you can use/implement your own loop.\ :raw-html-m2r:`<br/>`
An example is below:

.. code-block:: cpp

    #include "pressio/type_traits.hpp"
    #include "pressio/ode_solvers_nonlinear.hpp"
    #include "pressio/ode_advancers.hpp"
    #include "pressio/ode_steppers_implicit.hpp"

    int main()
    {
      MyRhsJacEvalClass mySystem(/*whatever args*/);

      namespace pode = pressio::ode;
      const auto scheme = pode::StepScheme::BDF1;
      auto stepper = pode::create_implicit_stepper(scheme, mySystem);

      // create a solver, here for simplicity we show the case where
      // for the types used, we can leverage pressio solvers
      using jacobian_t = typename problem_t::jacobian_type;
      using lin_solver_t = pressio::linearsolvers::Solver</*whatever*/, jacobian_t>;
      lin_solver_t linSolverObj;
      auto nonLinSolver = nonlinearsolvers::create_newton_raphson(stepper, linSolverObj);

      auto myState = mySystem.createState();
      // initialize myState to initial condition

      using time_type = typename decltype(mySystem)::independent_variable_type;
      const time_type t0 = 0.;
      const time_type dt = 0.1;
      const int num_steps = 22;
      pode::advance_n_steps(stepper, stateObj, t0, dt, num_steps, nonLinearSolver);
    }
