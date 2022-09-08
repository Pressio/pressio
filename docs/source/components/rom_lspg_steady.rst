
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
  /*impl defined*/ create_steady_problem(const TrialSpaceType & trialSpace,         (1)
                                         const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianActionMaskerType
    >
  /*impl defined*/ create_steady_problem(const TrialSpaceType & trialSpace,         (2)
					 const FomSystemType & fomSystem,
					 const ResidualMaskerType & rMasker,
					 const JacobianActionMaskerType & jaMasker);

  }}} // end namespace pressio::rom::lspg

- 1: overload for default OR hyper-reduced problem

- 2: overload for masked problem

Parameters
----------

* ``trialSpace``: trial subspace approximating the full space

* ``fomSystem``: full-order model instance

* ``rMasker``: masking operator to apply to the FOM residual

* ``jaMasker``: masking operator to apply to the result of the FOM jacobian action

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must model the ``TrialColumnSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialColumnSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType`` must model the ``SteadyFomWithJacobianAction`` `concept <rom_concepts/c6.html>`__

- ``ResidualMaskerType`` must model the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

- ``JacobianActionMaskerType`` must model the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

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

- The return value is an instance of a implementation-defined class
  representing a LSPG steady problem.
  This problem class is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class SteadyLspgProblemExpositionOnly
    {
      public:
        using state_type                = /* same as reduced_state in trialSpace*/;
        using residual_type             = /* impl defined*/;
        using jacobian_type	        = /* impl defined*/;

	state_type    createState() const;
        residual_type createResidual() const;
        jacobian_type createJacobian() const;
        void residualAndJacobian(const state_type & odeState,
	                         residual_type & R,
				 jacobian_type & J,
				 bool computeJacobian) const;
    };

.. important::

   Any steady LSPG problem satisfies the ``SystemWithFusedResidualAndJacobian``
   concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.

- the problem object will hold const-qualified references to the arguments
  ``trialSpace``, ``fomSystem``, ``rMasker``, ``jaMasker``, therefore
  NO copy of these objects occurs.

- All internal memory allocation needed for the implementation is
  performed inside the constructor of problem.


Using the problem
-----------------

The problem class satisfies the ``SystemWithFusedResidualAndJacobian`` concept
discussed `here <nonlinearsolvers_concepts/c2.html>`__.

To solve the problem, ... :red:`finish`
