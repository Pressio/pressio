
.. include:: ../mydefs.rst

Galerkin: Steady
================

Header: ``<pressio/rom_galerkin_steady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace galerkin{

  template<
    class TrialSpaceType,
    class FomSystemType>
  /*impl defined*/ create_steady_problem(const TrialSpaceType & trialSpace,         (1)
                                         const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType
    >
  /*impl defined*/ create_steady_problem(const TrialSpaceType & trialSpace,         (2)
					 const FomSystemType & fomSystem,
					 const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianActionMaskerType,
    class HyperReductionOperatorType
    >
  /*impl defined*/ create_steady_problem(const TrialSpaceType & trialSpace,         (3)
					 const FomSystemType & fomSystem,
					 const ResidualMaskerType & rMasker,
					 const JacobianActionMaskerType & jaMasker,
					 const HyperReductionOperatorType & hrOp);

  }}} // end namespace pressio::rom::galerkin

- 1: overload for default problem

- 2: overload for hyper-reduced problem

- 3: overload for masked problem

Parameters
----------

* ``trialSpace``: linear trial subspace approximating the FOM state space

* ``fomSystem``: full-order model instance

* ``hrOp``: hyper-reduction operator to apply to residual and jacobian action

* ``rMasker``: masking operator to apply to the FOM residual

* ``jaMasker``: masking operator to apply to the result of the FOM jacobian action

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must meet the ``TrialSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType`` must meet the ``SteadyFomWithJacobianAction`` `concept <rom_concepts/c6.html>`__

- ``HyperReductionOperatorType`` must meet the ``SteadyGalerkinHyperReductionOperator`` `concept <rom_concepts/c4.html>`__

- ``ResidualMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

- ``JacobianActionMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- all arguments passed to the function must be lvalues with a lifetime
  *longer* that that of the instantiated problem, i.e., they must be
  destructed *after* the problem goes out of scope

Mandates
~~~~~~~~

- the type representing the FOM state declared inside the ``TrialSpaceType``
  must be equal to that declared inside the ``FomSystemType`` class,
  i.e.: ``std::is_same<typename TrialSpaceType::full_state_type,
  typename FomSystemType::state_type >::value == true``

Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The functions construct and return an instance of an implementation-defined class
  representing a Galerkin steady problem. This problem class is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class SteadyGalerkinProblemExpositionOnly
    {
      public:
        using state_type                = /*same as reduced_state in TrialSpaceType*/;
        using residual_type             = /*impl defined*/;
        using jacobian_type	        = /*impl defined*/;

	state_type    createState() const;
        residual_type createResidual() const;
        jacobian_type createJacobian() const;
        void residualAndJacobian(const state_type & odeState,
	                         residual_type & R,
				 jacobian_type & J,
				 bool computeJacobian) const;
    };

.. important::

   Any steady Galerkin problem satisfies the ``SystemWithFusedResidualAndJacobian``
   concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.

- the problem object will hold const-qualified references to the arguments
  ``trialSpace``, ``fomSystem``, ``hrOp``, ``rMasker``, ``jaMasker``, therefore
  NO copy of these objects occurs.

- All internal memory allocation needed for the implementation is
  performed inside the constructor of problem.


Using the problem
-----------------

The problem class satisfies the ``SystemWithFusedResidualAndJacobian`` concept
discussed `here <nonlinearsolvers_concepts/c2.html>`__.

To solve the problem, you can use the Newton-Raphson solver from the pressio/nonlinear_solvers.
Or you can use/implement your own.

:red:`finish`
