
.. include:: ../mydefs.rst

Galerkin: Unsteady (implicit)
=============================

Header: ``<pressio/rom_galerkin_unsteady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace galerkin{

  template<
    class TrialSubspaceType,
    class FomSystemType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteadyexplicit::ComposableIntoDefaultProblem<
                TrialSubspaceType, FomSystemType>
  #endif
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (1)
						    const TrialSubspaceType & trialSubspace,
						    const FomSystemType & fomSystem);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class HyperReducerType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteadyexplicit::ComposableIntoHyperReducedProblem<
                TrialSubspaceType, FomSystemType, HyperReducerType>
  #endif
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (2)
						    const TrialSubspaceType & trialSubspace,
						    const FomSystemType & fomSystem,
						    const HyperReducerType & hyperReducer);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class MaskerType,
    class HyperReducerType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires unsteadyexplicit::ComposableIntoHyperReducedMaskedProblem<
                TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>
  #endif
  /*impl defined*/ create_unsteady_implicit_problem(ode::StepScheme schemeName,       (3)
						    const TrialSubspaceType & trialSubspace,
						    const FomSystemType & fomSystem,
						    const MaskerType & masker,
						    const HyperReducerType & hyperReducer);

  }}} // end namespace pressio::rom::galerkin

Description
~~~~~~~~~~~

Overload set to instantiate a default (1), hyper-reduced (2) or masked (3) problem
with *implicit* time integration.

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``schemeName``
     - enum value to choose the desired time integration scheme

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

- `rom::galerkin::unsteadyimplicit::ComposableIntoDefaultProblem <rom_concepts_implicit_galerkin/default.html>`__

- `rom::galerkin::unsteadyimplicit::ComposableIntoHyperReducedProblem <rom_concepts_implicit_galerkin/hr.html>`__

- `rom::galerkin::unsteadyimplicit::ComposableIntoHyperReducedMaskedProblem <rom_concepts_implicit_galerkin/masked.html>`__


Preconditions
~~~~~~~~~~~~~

.. _implicitGalerkinPreconditions:

1. ``schemeName`` must be an implicit scheme,
   see `this page <ode_steppers_implicit.html>`__ for the choices

2. all arguments passed to ``create_unsteady_implicit_problem`` must have a
   lifetime *longer* that that of the instantiated problem, i.e., they must be
   destructed *after* the problem instantiated goes out of scope

3. the trial subspace must be an admissible approximation
   of the specific full state/problem represented by the ``fomSystem`` instance


Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the return is an instance of class representing an *implicit* Galerkin problem.

    The return type is implementation defined, but guaranteed to
    model the ``SteppableWithAuxiliaryArgs`` concept discussed `here <ode_concepts/c7.html>`__.

- Purely syntactically, the problem class API is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      template<class TrialSubspaceType, ...> // exposition only
      class UnsteadyImplicitGalerkinProblemExpositionOnly
      {
	public:
          using state_type = typename TrialSubspaceType::reduced_state_type;
	  using independent_variable_type = /*same as declared inside your FomSystemType*/;

	  template<SolverType>
	  void operator()(state_type & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<independent_variable_type> & /**/,
			  SolverType & /**/);
      };


- any necessary memory allocation needed for the implementation
  occurs when the constructor of the class is called. However, we
  guarantee (for now) that the implementation only uses via *const references*
  (as opposed to copying) the arguments passed to ``create_unsteady_implicit_problem``.
  This is why it is critical to ensure :ref:`precondition 2 <implicitGalerkinPreconditions>`
  is satisfied.

Solve the problem
-----------------

An unsteady implicit Galerkin problem satisfies the "steppable with args"
concept discussed `here <ode_concepts/c7.html>`__.
Solving these problems can thus be done via the "advancers"
in pressio/ode to step forward the problem.

:red:`finish`
