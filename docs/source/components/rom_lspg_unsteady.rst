
.. include:: ../mydefs.rst

LSPG: Unsteady
==============

Header: ``<pressio/rom_lspg_unsteady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace lspg{

  template<
    class TrialSpaceType,
    class FomSystemType>
  auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    (1)
                               const TrialSpaceType & trialSpace,
                               const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    (2)
			       const TrialSpaceType & trialSpace,
			       const FomSystemType & fomSystem,
			       const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianMaskerType>
  auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    (3)
			       const TrialSpaceType & trialSpace,
			       const FomSystemType & fomSystem,
			       const ResidualMaskerType & rMasker,
			       const JacobianMaskerType & jMasker);

  }}} // end namespace pressio::rom::lspg

- 1: overload for default problem

- 2: overload for hyper-reduced problem

- 3: overload for masked problem

Parameters
----------

* ``schemeName``: enum value to set the desired *implicit* scheme to use

* ``trialSpace``: the linear trial subspace to approximate the full space

* ``fomSystem``: your full-order problem

* ``rMasker``: operator for masking the LSPG residual

* ``jMasker``: operator for masking the LSPG jacobian


Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must meet the ``TrialSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType`` must meet the ``SemiDiscreteFomWithJacobianAction`` `concept <rom_concepts/c2.html>`__.

- ``ResidualMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

- ``JacobianMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- ``schemeName`` must be one of ``pressio::ode::StepScheme::{BDF1, BDF2}``

- the trial space and system arguments passed to the function must be lvalues
  with a lifetime that is *longer* that that of the instantiated problem, i.e., they are
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

The returned problem object is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class UnsteadyLspgProblemExpositionOnly
    {
      public:
        // these nested type aliases must be here
        using independent_variable_type = /* same as in your system class */;
        using state_type                = /* same as your reduced_state type  */;
        using residual_type             = /* same as the right_hand_side_type in your system class */;
        using jacobian_type	        = /* same as in your basis_type */;

        template<class SolverType>
        void operator()(StateType & /**/,
			const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			pressio::ode::StepCount /**/,
			pressio::ode::StepSize<independent_variable_type> /**/,
			SolverType & /**/)

      !!!!!!!!!! missing create state !!!!!!!!!!!!

        residual_type createResidual() const;
        jacobian_type createJacobian() const;
        void residualAndJacobian(const state_type & odeState,
	                         residual_type & R,
				 jacobian_type & J,
				 bool computeJacobian) const;
    };


Using/solving the problem
-------------------------

The problem class  satisfies

- the ``SteppableWithAuxiliaryArgs`` concept discussed `here <ode_concepts/c7.html>`__.

- the ``SystemWithFusedResidualAndJacobian`` concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.

To solve the problem, you need to use a non-linear least squares solver
from the nonlinear_solvers and the "advance" functions to step forward.
Or you can use/implement your own loop.
An example is below:

:red:`finish`
