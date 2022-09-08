.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

Fully-discrete systems
======================

Header: ``<pressio/ode_steppers_implicit.hpp>``

Public namespace: ``pressio::ode``


Scope
-----

An "implicit stepper" in pressio is an abstraction that represents
"how" to take a step when applying an :orange:`implicit scheme`
to initial value problems expressable in the following *fully discrete* form

.. math::

    \boldsymbol{R}(\boldsymbol{y_{n}}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}, \qquad y(t_0) = y_0

Here, :math:`y_{n}` is the state at the n-th step,
:math:`t` is independent variable, and :math:`R` is the residual.

.. admonition:: :medium:`Recall the definition of implicit methods`

    :medium:`Implicit methods update the state by solving a (potentially nonlinear) system of equations involving the current, the predicted state and possibly previous states.`


API
---

.. code-block:: cpp

   namespace pressio{ namespace ode{

   template<std::size_t TotalNumberOfDesiredStates, class SystemType>
   /*impl defined*/ create_implicit_stepper(const SystemType & system);

   }} //end namespace pressio::ode


Parameters and templates
~~~~~~~~~~~~~~~~~~~~~~~~

* ``TotalNumberOfDesiredStates``: defines how many *totaly* states you need/want.
  For example, say that your arbitrary scheme needs :math:`y_n+1, y_n, y_n-1`.
  In such case, you use: ``TotalNumberOfDesiredStates = 3``.
  If you need three auxiliary states (beside) the main state to update,
  then use: ``TotalNumberOfDesiredStates = 4``.

* ``system``: an object that evaluates the discrete residual and jacobian

* ``massMatOperator``: the mass matrix operator


Constraints
~~~~~~~~~~~

* ``TotalNumberOfDesiredStates``: currently must be set to one of `{2, 3, 4}`

- ``SystemType`` must model the ``FullyDiscreteSystemWithJacobian`` `concept <ode_concepts/c5.html>`__.


Preconditions
~~~~~~~~~~~~~

- ``system`` must bind to an lvalue object whose lifetime is
  longer that that of the instantiated stepper, i.e., it is destructed
  *after* the stepper goes out of scope


Return value
~~~~~~~~~~~~

An instance of a pressio implicit stepper for doing this "arbitrary scheme".
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


:red:`describe here what is going on`

Use the stepper
---------------

The stepper created using the functions satisfies two concepts:

- the ``SteppableWithAuxiliaryArgs`` concept discussed `here <ode_concepts/c7.html>`__.

- the ``SystemWithFusedResidualAndJacobian`` concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.
