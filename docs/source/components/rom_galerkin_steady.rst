
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
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSpace,         (1)
                                         const FomSystemType & fomSystem);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class HyperReducerType>
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSpace,         (2)
					 const FomSystemType & fomSystem,
					 const HyperReducerType & hypReducer);

  template<
    class TrialSubspaceType,
    class FomSystemType,
    class ResidualMaskerType,
    class JacobianActionMaskerType,
    class HyperReducerType>
  /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSpace,         (3)
					 const FomSystemType & fomSystem,
					 const ResidualMaskerType & rMasker,
					 const JacobianActionMaskerType & jaMasker,
					 const HyperReducerType & hypReducer);

  }}} // end namespace pressio::rom::galerkin

Description
~~~~~~~~~~~

- 1: instantiates a default problem

- 2: instantiates a hyper-reduced problem

- 3: instantiates a masked problem

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``trialSpace``
     - trial subspace approximating the FOM state space

   * - ``fomSystem``
     - full-order model instance

   * - ``hypReducer``
     - hyper-reduction operator to apply to FOM residual and jacobian action

   * - ``rMasker``
     - masking operator to apply to the FOM residual

   * - ``jaMasker``
     - masking operator to apply to the result of the FOM jacobian action

Constraints
~~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``TrialSubspaceType``
     - must model the `PossiblyAffineTrialColumnSubspace concept <rom_concepts/c10.html>`__

   * - ``FomSystemType``
     - must model the `SteadyFomWithJacobianAction concept <rom_concepts/c6.html>`__

   * - ``HyperReducerType``
     - must model the `SteadyGalerkinHyperReducer concept <rom_concepts/c4.html>`__

   * - ``ResidualMaskerType``
     - must model the `TimeInvariantMasker concept <rom_concepts/c3.html>`__

   * - ``JacobianActionMaskerType``
     - must model the `TimeInvariantMasker concept <rom_concepts/c3.html>`__

Mandates
~~~~~~~~

- ``std::is_same<typename TrialSubspaceType::full_state_type,
  typename FomSystemType::state_type >::value == true``

Preconditions
~~~~~~~~~~~~~

- all arguments passed to the function must be lvalues with a lifetime
  *longer* that that of the instantiated problem, i.e., they must be
  destructed *after* the problem goes out of scope

Return value, Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- the overload set returns an instance of class representing a Galerkin steady problem

    The return type is implementation defined, but guaranteed
    to model the ``SystemWithFusedResidualAndJacobian``
    concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.


- Purely syntactically, the problem class API is:

  .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class SteadyGalerkinProblemSyntaxOnly
      {
	public:
	  using state_type    = /**/;
	  using residual_type = /**/;
	  using jacobian_type = /**/;

	  state_type    createState() const;
	  residual_type createResidual() const;
	  jacobian_type createJacobian() const;
	  void residualAndJacobian(const state_type & reducedState,
				   residual_type & R,
				   jacobian_type & J,
				   bool computeJacobian) const;
      };

  where:

  - the ``state_type`` typedef is same as ``typename TrialSubspaceType::reduced_state_type``,
    so it is either an Eigen vector or a Kokkos view (see `here <rom_concepts/c10.html>`__)

  - ``residual_type``, ``jacobian_type`` are the types of
    the *reduced* residual and jacobian:

    - if ``state_type`` is an Eigen vector, then ``residual_type``
      is an Eigen vector and ``jacobian_type`` is an Eigen dense matrix

    - if ``state_type`` is a Kokkos rank-1 view, then ``residual_type``
      is a Kokkos rank-1 view and ``jacobian_type`` is a Kokkos rank-2 view


- the instantiated problem object will only reference the arguments passed
  to ``trialSpace``, ``fomSystem``, ``hypReducer``, ``rMasker``, ``jaMasker``,
  therefore NO copy of these objects occurs.

- All internal memory allocation needed for the implementation is
  performed inside the constructor of problem.


Solve the problem
-----------------

Solving a steady Galerkin problem boils down to solving
a reduced dense system of (nonlinear) equations.
A Galerkin problem object exposes the operators
defining this reduced system. To solve it, one can use the
Newton-Raphson solver from the pressio/nonlinear_solvers,
or use/implement their own.

Using pressio nonlinear solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   int main()
   {
     namespace pls  = pressio::linearsolvers;
     namespace pnls = pressio::nonlinearsolvers;
     namespace pgal = pressio::rom::galerkin;

     // assuming trialSpace and fomSystem already exist
     auto problem = pgal::create_steady_problem(trialSpace, fomSystem);

     using jacobian_type   = typename decltype(problem)::jacobian_type;
     using linear_solver_t = pls::Solver<pls::iterative::LSCG, jacobian_type>;
     auto nonLinearSolver  = pnls::create_newton_raphson(problem, linear_solver_t{});

     auto reducedState = trialSpace.createReducedState();
     // set initial condition for reducedState somehow

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
     namespace pls  = pressio::linearsolvers;
     namespace pnls = pressio::nonlinearsolvers;
     namespace pgal = pressio::rom::galerkin;

     // assuming trialSpace and fomSystem already exist
     auto problem = pgal::create_steady_problem(trialSpace, fomSystem);

     CustomSolver mySolver;
     auto reducedState = trialSpace.createReducedState();
     // set initial condition for reducedState somehow
     mySolver.doSolve(problem, reducedState);
   }
