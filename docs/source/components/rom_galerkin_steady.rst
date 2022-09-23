
.. include:: ../mydefs.rst

Galerkin: Steady
================

Header: ``<pressio/rom_galerkin_steady.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{ namespace galerkin{

  template<
    class TrialSubspaceType,
    class FomSystemType>
  /*
    requires steady::ProjectableOnPossiblyAffineSubspace<FomSystemType, TrialSubspaceType>
  */
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (1)
                                         const FomSystemType & fomSystem);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class HyperReducerType>
  /*
    requires steady::HyperReduceableWith<FomSystemType, HyperReducerType, TrialSubspaceType>
  */
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (2)
					 const FomSystemType & fomSystem,
					 const HyperReducerType & hyperReducer);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class MaskerType,
    class HyperReducerType>
  /*
    requires steady::HyperReduceableAndMaskableWith<FomSystemType, MaskerType,
						    HyperReducerType, TrialSubspaceType>
  */
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (3)
					 const FomSystemType & fomSystem,
					 const maskerType & masker,
					 const HyperReducerType & hyperReducer);

  }}} // end namespace pressio::rom::galerkin

Description
~~~~~~~~~~~

Overload set to instantiate a default (1), hyper-reduced (2) or masked problem (3).

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

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

- `rom::galerkin::steady::ProjectableOnPossiblyAffineSubspace <rom_concepts_steady_galerkin/steady_gal_default.html>`__

- `rom::galerkin::steady::HyperReduceableWith <rom_concepts_steady_galerkin/steady_gal_hr.html>`__

- `rom::galerkin::steady::HyperReduceableAndMaskableWith <rom_concepts_steady_galerkin/steady_gal_hr_mask.html>`__


Mandates
~~~~~~~~

- ``std::is_same<typename TrialSubspaceType::full_state_type,
  typename FomSystemType::state_type >::value == true``

Preconditions
~~~~~~~~~~~~~

.. _steadyGalerkinPreconditions:

1. all arguments passed to ``create_steady_problem`` must have a
   lifetime *longer* that that of the instantiated problem, i.e., they must be
   destructed *after* the problem instantiated goes out of scope

2. the trial subspace must be an admissible approximation
   of the specific full state/problem represented by the ``fomSystem`` instance


Return value, Postconditions, Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the return is an instance of class representing a Galerkin steady problem

    The return type is implementation defined, but guaranteed
    to model the ``SystemWithFusedResidualAndJacobian``
    concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.


- Purely syntactically, the problem class API is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      template<class TrialSubspaceType, ...> // exposition only
      class SteadyGalerkinProblemSyntaxOnly
      {
        using galerkin_types = SteadyGalerkinTraits<
	  typename TrialSubspaceType::reduced_state_type>;  // exposition only

	public:
	  using state_type    = typename TrialSubspaceType::reduced_state_type;
	  using residual_type = typename galerkin_types::reduced_residual_type;
	  using jacobian_type = typename galerkin_types::reduced_jacobian_type;

	  state_type    createState() const;
	  residual_type createResidual() const;
	  jacobian_type createJacobian() const;
	  void residualAndJacobian(const state_type & /*reducedState*/,
				   residual_type    & /*reducedResidual*/,
				   jacobian_type    & /*reducedJacobian*/,
				   bool               /*computeJacobian*/) const;
      };

  where:

  - the ``state_type`` is an alias to the reduced state type of your trial subspace class

  - ``residual_type``, ``jacobian_type`` are the *reduced* residual and jacobian types
    which are defined based on what the reduced state type is. Specifically, we have:

    - if ``state_type`` is an Eigen vector, then ``residual_type``
      is an Eigen vector and ``jacobian_type`` is an Eigen dense matrix

    - if ``state_type`` is a Kokkos rank-1 view, then ``residual_type``
      is a Kokkos rank-1 view and ``jacobian_type`` is a Kokkos rank-2 view


- any necessary memory allocation needed for the implementation
  occurs when the constructor of the class is called. However, we
  guarantee (for now) that the implementation only uses via *const references*
  (as opposed to copying) the arguments passed to the ``create_steady_problem``.
  This is why it is critical to ensure :ref:`precondition 1 <steadyGalerkinPreconditions>`
  is satisfied.


Problem class: description
--------------------------

:red:`finish, we need to connect to this to the math pdf when we have it`


Solve the problem
-----------------

Solving a steady Galerkin problem practically means solving
a **reduced dense system of (nonlinear) equations**.
A Galerkin problem object created as shown above exposes
the operators that define such reduced system.
Since the reduced system, by definition, constitutes a *determined system of equations*,
to solve it one can use a Newton-Raphson solver either
from the pressio/nonlinear_solvers, or use/implement their own.

Using pressio nonlinear solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   int main()
   {
     namespace pls  = pressio::linearsolvers;
     namespace pnls = pressio::nonlinearsolvers;
     namespace pgal = pressio::rom::galerkin;

     // assuming trialSubspace and fomSystem already exist
     auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

     using jacobian_type   = typename decltype(problem)::jacobian_type;
     using linear_solver_t = pls::Solver<pls::iterative::LSCG, jacobian_type>;
     auto nonLinearSolver  = pnls::create_newton_raphson(problem, linear_solver_t{});

     auto reducedState = problem.createState();
     // set initial condition for reducedState somehow
     // set other parameters for the solver if needed
     nonLinearSolver.solve(problem, reducedState);
   }

Use your own solver
~~~~~~~~~~~~~~~~~~~

If you don't want to use the pressio solvers,
you can easily set up your own because the problem object
fully identify the system to solve.

.. code-block:: cpp

   class CustomSolver{
     // constructor as needed

     template<class T>
     void doSolve(const T & problem,
                  typename T::state_type & state)
     {
       auto R = problem.createResidual();
       auto J = problem.createJacobian();

       for (...){
         system.residualAndJacobian(state, R, J, true);
         // do something, update state, etc.
       }
     }
   };

   int main()
   {
     namespace pgal = pressio::rom::galerkin;

     // assuming trialSubspace and fomSystem already exist
     auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

     CustomSolver mySolver;
     auto reducedState = problem.createState();
     // set initial condition for reducedState somehow
     mySolver.doSolve(problem, reducedState);
   }
