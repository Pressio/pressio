
.. include:: ../mydefs.rst

LSPG: Steady
============

Header: ``<pressio/rom_lspg_steady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace lspg{

  template<
    class TrialSpaceType,
    class FomSystemType>
  auto create_steady_problem(const TrialSpaceType & trialSpace,         (1)
                             const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianActionMaskerType
    >
  auto create_steady_problem(const TrialSpaceType & trialSpace,         (2)
			     const FomSystemType & fomSystem,
			     const ResidualMaskerType & rMasker,
			     const JacobianActionMaskerType & jaMasker);

  }}} // end namespace pressio::rom::lspg

- (1): overload for default OR hyper-reduced problem

- (2): overload for masked problem

Parameters
----------

* ``trialSpace``: the linear trial subspace to approximate the full space

* ``fomSystem``: your full-order problem

* ``rMasker``: operator for masking the FOM residual

* ``jaMasker``: operator for masking the result of the FOM jacobian action

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must meet the ``TrialColumnSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialColumnSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType`` must meet the ``SteadyFomWithJacobianAction`` `concept <rom_concepts/c6.html>`__

- ``ResidualMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

- ``JacobianActionMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- all arguments passed to the function must be lvalues with a lifetime
  that is *longer* that that of the instantiated problem, i.e., they are
  destructed *after* the problem goes out of scope

- the ``trialSpace`` must represent a space compatible with the ``fomSystem``

Mandates
~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

An instance of a implementation-defined class that represents a Galerkin steady problem.

Postconditions and side effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The returned problem object exposes this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class SteadyLspgProblemExpositionOnly
    {
      public:
        using state_type                = /* same as reduced_state in trialSpace*/;
        using residual_type             = /* impl defined*/;
        using jacobian_type	        = /* impl defined*/;

        residual_type createResidual() const;
        jacobian_type createJacobian() const;
        void residualAndJacobian(const state_type & odeState,
	                         residual_type & R,
				 jacobian_type & J,
				 bool computeJacobian) const;
    };

Use the problem
---------------

The problem class satisfies the ``SystemWithFusedResidualAndJacobian`` concept
discussed `here <nonlinearsolvers_concepts/c2.html>`__.

To solve the problem, ... :red:`finish`
