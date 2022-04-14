.. role:: raw-html-m2r(raw)
   :format: html

implicit steppers (discrete-time systems)
=========================================

.. note::

    Defined in: ``<pressio/ode_steppers_implicit.hpp>``

    Public namespace: ``pressio::ode``

Overview
--------

Provides functionalities to create steppers for implicit methods.
Recall that implicit methods update the state of a system
by solving a system of equations involving both the current and next state.
An implicit stepper is an object that knows how to do one such *implicit* step.

This page describes functionalities applicable to any system
expressible in a *discrete-time* form

.. math::

    \boldsymbol{R}(\boldsymbol{y_{n}}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}

Here, :math:`y_{n}` is the state at the n-th time step,
:math:`t` is time, and :math:`R` is the residual.


API, Parameters and Requirements
--------------------------------

\todo FINISH

.. code-block:: cpp

    template<int num_states, class StateType, class SystemType>
    auto create_arbitrary_stepper(const StateType & state,
                                  SystemType && system);

*
  ``num_states``: total number of states you need.

*
  ``state``:

  * self-explanatory

  * must be copy constructible and the following condition must be true:
    ``std::is_same<StateType, typename SystemType::state_type>::value == true``

*
  ``system``:

  * problem instance knowing how to create and compute the residual and Jacobian.

  * Must conform to:

    .. code-block:: cpp

        class ValidDiscreteTimeSystem
        {
        using scalar_type = /* ... */;
        using state_type  = /* ... */;
        using discrete_time_residual_type = /* ... */;
        using discrete_time_jacobian_type = /* ... */;

        discrete_time_residual_type createDiscreteTimeResidual() const;
        discrete_time_jacobian_type createDiscreteTimeJacobian() const;

        // overload accepting 1 auxiliary state
        template<class StepCountType>
        void discreteTimeResidual(StepCountType,
                                 scalar_type time,
                                 scalar_type dt,
                                 discrete_time_residual_type &,
                                 const state_type &) const;


        // overload accepting 2 auxiliary states
        template<class StepCountType>
        void discreteTimeResidual(StepCountType,
                                   scalar_type time,
                                   scalar_type dt,
                                   discrete_time_residual_type &,
                                   const state_type &,
                                   const state_type &) const;


        // overload accepting 1 auxiliary state
        template<class StepCountType>
        void discreteTimeJacobian(StepCountType,
                                   scalar_type time,
                                   scalar_type dt,
                                   discrete_time_jacobian_type &,
                                   const state_type &) const;


        // overload accepting 2 auxiliary states
        template<class StepCountType>
        void discreteTimeJacobian(StepCountType,
                                   scalar_type time,
                                   scalar_type dt,
                                   discrete_time_jacobian_type &,
                                   const state_type &,
                                   const state_type &) const;

        };
