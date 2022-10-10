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
   #ifdef PRESSIO_ENABLE_CXX20
     requires SystemWithRhsAndJacobian<mpl::remove_cvref_t<SystemType>>
   #endif
   /*impl defined*/ create_implicit_stepper(StepScheme odeScheme,                   (1)
                                            SystemType && system);

   template<class SystemType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires SystemWithRhsJacobianMassMatrix<mpl::remove_cvref_t<SystemType>>
   #endif
   /*impl defined*/ create_implicit_stepper(StepScheme odeScheme,                   (2)
                                            SystemType && system);

   }} //end namespace pressio::ode


Description
~~~~~~~~~~~

Overload set to instantiate a stepper without (1) or with (2) mass matrix.

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``odeScheme``
     - the target stepping scheme

   * - ``system``
     - problem instance

Constraints
~~~~~~~~~~~

Each overload is associated with a set of constraints.
With C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:

- `pressio::ode::SystemWithRhsAndJacobian <ode_concepts_system/rhs_jac.html>`__

- `pressio::ode::SystemWithRhsJacobianMassMatrix <ode_concepts_system/rhs_jac_massmatrix.html>`__


Preconditions
~~~~~~~~~~~~~

- ``odeScheme`` must be one of ``pressio::ode::StepScheme::{BDF1, BDF2}``.

- if ``system`` does *not* bind to a temporary object,
  it must bind to an lvalue object whose lifetime is
  *longer* that that of the instantiated stepper, i.e., it is destructed
  *after* the stepper goes out of scope


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

- if you pass a an rvalue "problem" object, the constructor of the stepper
  will try to use move semantics. If move semantics are implemented, the temporary
  is moved from and no new memory allocation occurs. If move semantics fallback to copying,
  then a copy of the original object is made. Either way, the stepper instantiated
  directly manages the lifetime of the object.

- if you pass a an lvalue "problem" object, the constructor of the stepper
  will create a *reference to the* object. In such case, the stepper
  object does NOT manage the lifetime of it. You are responsible to
  ensure that the ``system`` is destructed *after* the stepper goes out of scope

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
      MySystemClass mySystem(/*whatever args*/);

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
