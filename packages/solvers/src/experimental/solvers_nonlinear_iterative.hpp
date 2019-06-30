
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP

#include <type_traits>
// #include "containers_ConfigDefs.hpp"
#include "solvers_nonlinear_base.hpp"
#include "solvers_nonlinear_iterative_helper.hpp"
// #include "solvers_nonlinear_factory.hpp"


namespace rompp{
namespace solvers{


struct NonLinearSolvers; // Fwd declaration

/**
 * Implements a non linear solver bases on a linear iterative solver.
 */
template <
  typename PolicyT,
  typename LSolverT
>
class NonLinearIterativeSolver
  : public NonLinearIterativeSolverHelper,
    public NonLinearSolverBase<
      NonLinearIterativeSolver<
        PolicyT,
        LSolverT
      >
    >{

    friend NonLinearSolvers;
    typedef NonLinearSolverBase<
      NonLinearIterativeSolver<PolicyT, LSolverT>> base_type;

  public:

    /**
    * Implements the solve method for a non linear solver.
    */
    template <
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT
    >
    auto solve_(const SystemT& sys, const VectorT& b) 
    -> decltype( 
          PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b,
               1, 1, 1e-2, 1e-2) ) {

      double tolerance = this->getTolerance();
      double nonLinearTolerance = this->getNonLinearTolerance();

      containers::default_types::uint maxIterations = this->getMaxIterations();
      containers::default_types::uint maxNonLinearIterations = this->getMaxNonLinearIterations();

      return PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b,
        			 maxIterations, maxNonLinearIterations,
        			 tolerance, nonLinearTolerance);
    }
    //--------------------------------------------------------------

protected:

  NonLinearIterativeSolver()
    : NonLinearIterativeSolverHelper(), base_type() {}

};

} // end namespace solvers
}//end namespace rompp
#endif
