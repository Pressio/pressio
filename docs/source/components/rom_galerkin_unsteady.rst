
.. include:: ../mydefs.rst

Galerkin: Unsteady
==================

Header: ``<pressio/rom_galerkin_unsteady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace galerkin{

  //
  // overload set for explicit in time
  //
  template<
    class TrialSpaceType,
    class FomSystemType>
  auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,    (1)
					const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,    (2)
					const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class RhsMaskerType,
    class HyperReductionOperatorType>
  auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,    (3)
					const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const RhsMaskerType & rhsMasker,
					const HyperReductionOperatorType & hrOp);

  //
  // overload set for implicit in time
  //
  template<
    class TrialSpaceType,
    class FomSystemType>
  auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,    (4)
					const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,    (5)
					const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const HyperReductionOperatorType & hrOp);

  }}} // end namespace pressio::rom::galerkin

- 1,2,3: overload for (1) default, (2) hyper-reduced and (3) masked problem with *explicit* time integration

- 4,5,6: overload for (4) default, (5) hyper-reduced and (6) masked problem with *implicit* time integration

Parameters
----------

* ``schemeName``: enum value to set the desired time integration scheme to use

* ``trialSpace``: the linear trial subspace to approximate the full space

* ``fomSystem``: your full-order problem

* ``hrOp``: operator to left-multiply the hyper-reduced FOM rhs

* ``rhsMasker``: operator for masking the FOM right hand side

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must meet the ``TrialSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType``:

  - for 1,2,3, it must meet the ``SemiDiscreteFom`` `concept <rom_concepts/c1.html>`__.

  - for 4,5,6, it must meet the ``SemiDiscreteFomWithJacobianAction`` `concept <rom_concepts/c2.html>`__.

- ``HyperReductionOperatorType``:

  - for 2,3, it must meet the ``UnsteadyGalerkinRhsHyperReductionOperator`` `concept <rom_concepts/c4b.html>`__

  - for 5,6, it must meet the ``UnsteadyGalerkinRhsAndJacobianHyperReductionOperator`` `concept <rom_concepts/c4b.html>`__

- ``RhsMaskerType`` must meet the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- for 1,2,3: ``schemeName`` must be an explicit scheme,
  see `this page <ode_steppers_explicit.html>`__ for the choices

- for 4,5,6: ``schemeName`` must be an implicit scheme,
  see `this page <ode_steppers_implicit.html>`__ for the choices

- the trial space, system, masker, and hyper-reduction operator arguments passed
  to the function must be lvalues with a lifetime *longer* that that of
  the instantiated problem, i.e., they are destructed *after* the problem goes out of scope

- the ``trialSpace`` must represent a space compatible with the ``fomSystem``

Mandates
~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

An instance of a implementation-defined class that represents a Galerkin unsteady explicit problem.

Postconditions and side effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you choose an explicit problem, the returned problem object is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class UnsteadyExplicitGalerkinProblemExpositionOnly
    {
      public:
	using state_type                = /* same as the reduced_state_type */;
	using independent_variable_type = /* same as your FOM time_type */;

	void operator()(StateType & /**/,
			const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			pressio::ode::StepCount /**/,
			const pressio::ode::StepSize<IndVarType> & /**/);
    };

If you choose an implicit problem, the returned problem object is guaranteed to expose this API:

.. code-block:: cpp

    // This is not the actual class, it just describes the API
    class UnsteadyExplicitGalerkinProblemExpositionOnly
    {
      public:
	using state_type                = /* same as the reduced_state_type */;
	using independent_variable_type = /* same as your FOM time_type */;

	template<SolverType>
	void operator()(StateType & /**/,
			const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			pressio::ode::StepCount /**/,
			const pressio::ode::StepSize<IndVarType> & /**/,
			SolverType & /**/);
    };


Using/solving the problem
-------------------------

Two key things to notice here are:

- an unsteady explicit Galerkin problem satisfies the "steppable" concept
  discussed `here <ode_concepts/c6.html>`__

- an unsteady implicit Galerkin problem satisfies the "steppable with args" concept
  discussed `here <ode_concepts/c7.html>`__

so one can use the "advancers" in pressio/ode to step forward the problem.
Example below:

:red:`finish`
