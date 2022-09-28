
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
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteady::ComposableIntoDefaultProblem<
               TrialSubspaceType, FomSystemType>
  #endif
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,         (1)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteady::ComposableIntoHyperReducedProblem<
               TrialSubspaceType, FomSystemType, HyperReductionOperatorType>
  #endif
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,          (2)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem,
					   const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class MaskerType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteady::ComposableIntoMaskedProblem<
               TrialSubspaceType, FomSystemType, MaskerType>
  #endif
  /*impl defined*/ create_unsteady_problem(ode::StepScheme schemeName,          (3)
					   const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem,
					   const MaskerType & masker);

  template<
    std::size_t TotalNumberOfStencilStates,
    class TrialSpaceType,
    class FomSystemType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires FullyDiscreteSystemWithJacobianAction<
                FomSystemType, TotalNumberOfDesiredStates, TrialSubspaceType>
  #endif
  /*impl defined*/ create_unsteady_problem(const TrialSpaceType & trialSpace,   (4)
					   const FomSystemType & fomSystem);

  }}} // end namespace pressio::rom::lspg

Description
~~~~~~~~~~~

Overload set to instantiate a default (1), hyper-reduced (2), masked (3) or "user-defined" (4) problem.

Non-type Template Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``TotalNumberOfStencilStates``: total number of desired states needed
  to define your scheme. Applicable only to overload 4.

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``schemeName``
     - enum value to set the desired *implicit* scheme to use

   * - ``trialSubspace``
     - trial subspace approximating the FOM state space

   * - ``fomSystem``
     - full-order model instance

   * - ``hyperReducer``
     - hyper-reduction operator

   * - ``masker``
     - masking operator

Constraints
~~~~~~~~~~~

Each overload is associated with a set of constraints.
If we could use C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet use C++20, the constraints are
currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:

- `rom::lspg::unsteady::DefaultProblem <rom_concepts_unsteady_lspg/default.html>`__

- `rom::lspg::unsteady::HyperReducedProblem <rom_concepts_unsteady_lspg/hyperreduced.html>`__

- `rom::lspg::unsteady::MaskedProblem <rom_concepts_unsteady_lspg/masked.html>`__

- `rom::FullyDiscreteSystemWithJacobianAction <rom_concepts_foms/fully_discrete_with_jac_action.html>`__

Mandates
~~~~~~~~

- ``std::is_same<typename TrialSubspaceType::full_state_type,
  typename FomSystemType::state_type >::value == true``

Preconditions
~~~~~~~~~~~~~

.. _unsteadyLspgPreconditions:

1. ``schemeName`` must be one of ``pressio::ode::StepScheme::{BDF1, BDF2}``

2. all arguments passed to ``create_steady_problem`` must have a
   lifetime *longer* that that of the instantiated problem, i.e., they must be
   destructed *after* the problem instantiated goes out of scope

2. the trial subspace must be an admissible approximation
   of the specific full state/problem represented by the ``fomSystem`` instance

Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the return value is an instance of class representing an unsteady LSPG problem.

    The return type is implementation defined, but guaranteed to
    model the ``SteppableWithAuxiliaryArgs`` concept
    discussed `here <ode_concepts/c7.html>`__.

  This means that the purely syntactical API of the problem class is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class UnsteadyLspgProblemExpositionOnly
      {
	public:
	  using independent_variable_type = /*see description below*/;
	  using state_type                = /*see description below*/;

	  template<class SolverType>
	  void operator()(state_type & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  pressio::ode::StepSize<independent_variable_type> /**/,
			  SolverType & /**/)
      };

  where:

  - ``state_type`` aliases the reduced state type of your trial subspace class

  - ``independent_variable_type`` same as the nested type alias in ``FomSystemType``

- any necessary memory allocation needed for the implementation
  occurs when the constructor of the class is called. However, we
  guarantee (for now) that the implementation only uses via *const references*
  (as opposed to copying) the arguments passed to the ``create_steady_problem``.
  This is why it is critical to ensure :ref:`precondition 2 <unsteadyLspgPreconditions>`
  is satisfied.

Solve the problem
-------------------------

To solve the problem, you need a non-linear least squares solver
from the nonlinear_solvers and the "advance" functions to step forward.
Or you can use/implement your own loop.
A representative snippet is below:

:red:`finish`
