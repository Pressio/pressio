
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
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,       (1)
						    const TrialSpaceType & trialSpace,
						    const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,       (2)
						    const TrialSpaceType & trialSpace,
						    const FomSystemType & fomSystem,
						    const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class RhsMaskerType,
    class HyperReductionOperatorType>
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,       (3)
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
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (4)
						    const TrialSpaceType & trialSpace,
						    const FomSystemType & fomSystem);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class HyperReductionOperatorType>
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (5)
						    const TrialSpaceType & trialSpace,
						    const FomSystemType & fomSystem,
						    const HyperReductionOperatorType & hrOp);

  template<
    class TrialSpaceType,
    class FomSystemType,
    class RhsMaskerType,
    class JacobianActionMaskerType,
    class HyperReductionOperatorType>
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (6)
						    const TrialSpaceType & trialSpace,
						    const FomSystemType & fomSystem,
						    const RhsMaskerType & rhsMasker,
						    const JacobianActionMaskerType & jaMasker,
						    const HyperReductionOperatorType & hrOp);

  }}} // end namespace pressio::rom::galerkin

- 1,2,3: overloads for (1) default, (2) hyper-reduced and (3) masked problem with *explicit* time integration

- 4,5,6: overloads for (4) default, (5) hyper-reduced and (6) masked problem with *implicit* time integration

Parameters
~~~~~~~~~~

* ``schemeName``: enum value to choose the desired time integration scheme

* ``trialSpace``: trial subspace approximating the FOM state space

* ``fomSystem``: full-order model instance

* ``hrOp``: hyper-reduction operator

* ``rhsMasker``: masking operator to apply to the FOM right hand side

* ``jaMasker``: masking operator to apply to the result of the FOM jacobian action

Constraints
~~~~~~~~~~~

- ``TrialSpaceType`` must model the ``TrialColumnSubspace`` `concept <rom_concepts/c7.html>`__
  or ``AffineTrialColumnSubspace`` `concept <rom_concepts/c8.html>`__

- ``FomSystemType``:

  - for 1,2,3: it must model the ``SemiDiscreteFom`` `concept <rom_concepts/c1.html>`__.

  - for 4,5,6: it must model the ``SemiDiscreteFomWithJacobianAction`` `concept <rom_concepts/c2.html>`__.

- ``HyperReductionOperatorType``:

  - for 2,3: it must model the ``ExplicitGalerkinHyperReducer`` `concept <rom_concepts/c4b.html>`__

  - for 5,6: it must model the ``ImplicitGalerkinHyperReducer`` `concept <rom_concepts/c4c.html>`__

- ``RhsMaskerType`` and ``JacobianActionMaskerType`` must both meet
  the ``TimeInvariantMasker`` `concept <rom_concepts/c3.html>`__

Preconditions
~~~~~~~~~~~~~

- 1,2,3: ``schemeName`` must be an explicit scheme,
  see `this page <ode_steppers_explicit.html>`__ for the choices

- 4,5,6: ``schemeName`` must be an implicit scheme,
  see `this page <ode_steppers_implicit.html>`__ for the choices

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

  - for 3,6: ``std::is_same<
    typename RhsMaskerType::operand_type,
    typename FomSystemType::right_hand_side_type>::value == true``

  - for 6: let ``fom_jac_action_result_type`` the type of the
    result of applying the FOM Jacobian to the basis, then the following must hold:
    ``std::is_same<typename JacobianActionMaskerType::operand_type,
    fom_jac_action_result_type>::value == true``

- the hyper-reduction operato must "be compatible" with
  the masker, so we must have:

  - for 3,6: ``std::is_same<
    typename HyperReductionOperatorType::right_hand_side_operand_type,
    typename RhsMaskerType::result_type>::value == true``

  - for 6: the following must hold: ``std::is_same<
    typename HyperReductionOperatorType::jacobian_action_operand_type,
    typename JacobianActionMaskerType::result_type>::value == true``


Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the overload set 1,2,3 returns an instance of
  a class representing an *explicit* Galerkin problem.

    The return type is implementation defined, but guaranteed to
    model the ``Steppable`` concept discussed `here <ode_concepts/c6.html>`__.

  This means that the purely syntactical API of the problem class is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class UnsteadyExplicitGalerkinProblemExpositionOnly
      {
	public:
	  using state_type                = /*same as the reduced_state_type in TrialSpaceType*/;
	  using independent_variable_type = /*same as declared inside your FomSystemType*/;

	  void operator()(state_type & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<independent_variable_type> & /**/);
      };


- the overload set 4,5,6 returns an instance of
  a class representing an *implicit* Galerkin problem.

    The return type is implementation defined, but guaranteed to
    model the ``SteppableWithAuxiliaryArgs`` concept discussed `here <ode_concepts/c7.html>`__.

  This means that the purely syntactical API of the problem class is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class UnsteadyImplicitGalerkinProblemExpositionOnly
      {
	public:
	  using state_type                = /*same as the reduced_state_type in TrialSpaceType*/;
	  using independent_variable_type = /*same as declared inside your FomSystemType*/;

	  template<SolverType>
	  void operator()(state_type & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<independent_variable_type> & /**/,
			  SolverType & /**/);
      };


- for all overloads, the problem object will hold const-qualified references
  to the arguments ``trialSpace``, ``fomSystem``, ``hrOp``, ``rhsMasker``, ``jaMasker``,
  therefore NO copy of these objects occurs.

- All internal memory allocation needed for the implementation is
  performed inside the constructor of problem.


Using/solving the problem
-------------------------

Two key things to notice here are:

- an unsteady explicit Galerkin problem satisfies the "steppable" concept
  discussed `here <ode_concepts/c6.html>`__

- an unsteady implicit Galerkin problem satisfies the "steppable with args" concept
  discussed `here <ode_concepts/c7.html>`__

Solving these problems can thus be done via the "advancers"
in pressio/ode to step forward the problem.

:red:`finish`
