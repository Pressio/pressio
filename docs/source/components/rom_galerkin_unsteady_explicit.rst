
.. include:: ../mydefs.rst

Galerkin: Unsteady (explicit)
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
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,         (1)
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
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,         (2)
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
  /*impl defined*/ create_unsteady_explicit_problem(ode::StepScheme schemeName,         (3)
						    const TrialSubspaceType & trialSubspace,
						    const FomSystemType & fomSystem,
						    const MaskerType & masker,
						    const HyperReducerType & hyperReducer);

  }}} // end namespace pressio::rom::galerkin

Description
~~~~~~~~~~~

Overload set to instantiate a default (1), hyper-reduced (2) or masked (3) problem
with *explicit* time integration.

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
With C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:


- `rom::galerkin::unsteadyexplicit::ComposableIntoDefaultProblem <rom_concepts_explicit_galerkin/default.html>`__

- `rom::galerkin::unsteadyexplicit::ComposableIntoHyperReducedProblem <rom_concepts_explicit_galerkin/hr.html>`__

- `rom::galerkin::unsteadyexplicit::ComposableIntoHyperReducedMaskedProblem <rom_concepts_explicit_galerkin/masked.html>`__

Preconditions
~~~~~~~~~~~~~

.. _explicitGalerkinPreconditions:

1. ``schemeName`` must be an explicit scheme,
   see `this page <ode_steppers_explicit.html>`__ for the choices

2. all arguments passed to ``create_unsteady_explicit_problem`` must have a
   lifetime *longer* that that of the instantiated problem, i.e., they must be
   destructed *after* the problem instantiated goes out of scope

3. the trial subspace must be an admissible approximation
   of the specific full state/problem represented by the ``fomSystem`` instance


Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the return is an instance of class representing an *explicit* Galerkin problem.

    The return type is implementation defined, but guaranteed to
    model the ``Steppable`` concept discussed `here <ode_concepts/c6.html>`__.


- Purely syntactically, the problem class API is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      template<class TrialSubspaceType, ...> // exposition only
      class UnsteadyExplicitGalerkinProblemExpositionOnly
      {
	public:
          using state_type = typename TrialSubspaceType::reduced_state_type;
	  using independent_variable_type = /*same as in your FomSystemType*/;

	  void operator()(state_type & /**/,
			  const pressio::ode::StepStartAt<independent_variable_type> & /**/,
			  pressio::ode::StepCount /**/,
			  const pressio::ode::StepSize<independent_variable_type> & /**/);
      };


- any necessary memory allocation needed for the implementation
  occurs when the constructor of the class is called. However, we
  guarantee (for now) that the implementation only uses via *const references*
  (as opposed to copying) the arguments passed to ``create_unsteady_explicit_problem``.
  This is why it is critical to ensure :ref:`precondition 2 <explicitGalerkinPreconditions>`
  is satisfied.


Solve the problem
-----------------

An unsteady explicit Galerkin problem satisfies the "steppable" concept
discussed `here <ode_concepts/c6.html>`__, therefore you can use
the step it to move forward (or backward) in time.
To do so, you can use the advancers" in ``pressio/ode`` as shown below,
or you can set up your own time stepping loop.

.. code-block:: cpp

   #include "pressio/rom_subspaces.hpp"
   #include "pressio/rom_galerkin_unsteady.hpp"

   template<class ReducedStateType>
   struct ReducedStateObserver
   {
     template<class TimeType>
     void operator()(pressio::ode::StepCount stepIn,
		     TimeType currentTime,
		     const ReducedStateType & romState) const
     {
       // you are given the step count, the time
       // and the corresponding rom state
       // so you can observe, store, etc as needed
     }
   };

   int main()
   {
     namespace pode = pressio::ode;
     namespace pgal = pressio::rom::galerkin;

     /*assuming:
      - trialSubspace : created somehow
      - fomSystem     : an instance of the full-order model
     */
     auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

     const auto odeScheme = pode::StepScheme::RungeKutta4;
     auto problem = pgal::create_unsteady_explicit_problem(odeScheme, trialSubpace, fomSystem);

     auto romState = trialSubspace.createReducedState();
     // set reduced state initial condition somehow

     using time_type = typename fom_t::time_type;
     const time_type dt = /*set time step size*/;
     ReducedStateObserver<decltype(romState)> observer;
     pode::advance_n_steps(problem,              /*the problem is our steppable object*/
                           romState,             /*the state to evolve in time*/
			   time_type{0},         /*start at time = 0*/
			   dt,                   /*the time step size*/
			   pode::StepCount(100), /*how many steps to take*/
			   observer);            /*an observer to monitor the solution*/
   }
