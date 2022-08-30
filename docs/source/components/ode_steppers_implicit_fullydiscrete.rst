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

..
   Usage
   -----

   The main steps are:

   1. you define your problem
   2. you create a pressio stepper
   3. you use the stepper to advance

   1. Define your problem
   ----------------------

   Your problem class must satisfy the ``FullyDiscreteSystemWithJacobian`` `concept <ode_concepts/c5.html>`__.

   Practical Considerations
   ^^^^^^^^^^^^^^^^^^^^^^^^

   :red:`describe here what is going on`

   2. Instantiate a stepper
   ------------------------

   .. code-block:: cpp

       template<int num_states, class SystemType>
       auto create_implicit_stepper(SystemType && system);

   * ``num_states``

     * this is the *TOTAL* number of states you need for your scheme.
       For example, say that your arbitrary scheme needs :math:`y_n+1, y_n, y_n-1`.
       In such case, you use: ``num_states = 3``.
       If you need three auxiliary states (beside) the main state to update,
       then use: ``num_states = 4``.

   * ``system``

     * problem instance satisfying the ``FullyDiscreteSystemWithJacobian`` `concept <ode_concepts/c5.html>`__.



   The stepper class API
   ^^^^^^^^^^^^^^^^^^^^^

   The returned stepper object has this API:

   .. code-block:: cpp

       // This is not the actual class, it just describes the API
       class Stepper
       {
	 public:
	   // these nested type aliases must be here
	   using independent_variable_type = /* same as in your system class */;
	   using state_type                = /* same as in your system class */;
	   using residual_type             = /* equal to discrete_residual_type in your system class */;
	   using jacobian_type	        = /* equal to discrete_jacobian_type in your system class */;

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

   |

   3. Use the stepper
   ------------------

   You use the stepper similarly to what shown for
   the `implicit steppers of a semi-discrete problems <ode_steppers_implicit_semidiscrete.html>`__.
