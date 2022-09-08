
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
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,       (1)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,        (2)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem,
					   const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianMaskerType>
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,        (3)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem,
					   const ResidualMaskerType & rMasker,
					   const JacobianMaskerType & jMasker);

  template<
    std::size_t TotalNumberOfDesiredStates,
    class TrialSpaceType,
    class FomSystemType>
  /*impl defined*/ create_unsteady_problem(const TrialSpaceType & trialSpace, (4)
					   const FomSystemType & fomSystem);

  }}} // end namespace pressio::rom::lspg

- 1: overload for default problem

- 2: overload for hyper-reduced problem

- 3: overload for masked problem

- 4: overload for an "arbitrary" problem

Templates and Parameters
------------------------

* ``schemeName``: enum value to set the desired *implicit* scheme to use

* ``trialSpace``: trial subspace to approximate the full space

* ``fomSystem``: full-order model instance

* ``rMasker``: operator for masking the LSPG residual

* ``jMasker``: operator for masking the LSPG jacobian

* ``TotalNumberOfDesiredStates``: total number of desired states needed to define your scheme

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must model the ``TrialColumnSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialColumnSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType``:

  - for 1,2,3: must model the ``SemiDiscreteFomWithJacobianAction`` `concept <rom_concepts/c2.html>`__.

  - for 4: must model the ``FullyDiscreteFomWithJacobianAction`` `concept <rom_concepts/c5.html>`__

- ``ResidualMaskerType`` must model the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

- ``JacobianMaskerType`` must model the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- ``schemeName`` must be one of ``pressio::ode::StepScheme::{BDF1, BDF2}``

- all arguments passed to the function must be lvalues with a lifetime
  *longer* that that of the instantiated problem, i.e., they must be
  destructed *after* the problem goes out of scope

Mandates
~~~~~~~~

- the type representing the FOM state declared inside the ``TrialSpaceType``
  must be equal to that declared inside the ``FomSystemType`` class,
  i.e.: ``std::is_same<typename TrialSpaceType::full_state_type,
  typename FomSystemType::state_type >::value == true``

- the masking operators must be compatible with the FOM types,
  so we must have:

  - ``std::is_same<
    typename ResidualMaskerType::operand_type,
    typename FomSystemType::right_hand_side_type>::value == true``

  - let ``fom_jac_action_result_type`` the type of the
    result of applying the FOM Jacobian to the basis, then the following must hold:
    ``std::is_same<typename JacobianActionMaskerType::operand_type,
    fom_jac_action_result_type>::value == true``

Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- An instance of a implementation-defined class representing a LSPG unsteady problem.
  This problem class is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class UnsteadyLspgProblemExpositionOnly
    {
      public:
        // these nested type aliases must be here
        using independent_variable_type = /* same as in your system class */;
        using state_type                = /* same as your reduced_state type  */;

        template<class SolverType>
        void operator()(StateType & /**/,
			const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			pressio::ode::StepCount /**/,
			pressio::ode::StepSize<independent_variable_type> /**/,
			SolverType & /**/)
    };

.. important::

   Any unsteady LSPG problem satisfies the ``SteppableWithAuxiliaryArgs``
   concept discussed `here <ode_concepts/c7.html>`__.

- for all overloads, the problem object will hold const-qualified references
  to the arguments ``trialSpace``, ``fomSystem``, ``hrOp``, ``rhsMasker``, ``jaMasker``,
  therefore NO copy of these objects occurs.

- All internal memory allocation needed for the implementation is
  performed inside the constructor of problem.

Using/solving the problem
-------------------------

To solve the problem, you need a non-linear least squares solver
from the nonlinear_solvers and the "advance" functions to step forward.
Or you can use/implement your own loop.
An example is below:

:red:`finish`
