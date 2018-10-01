
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP

#include <type_traits>
// #include "core_ConfigDefs.hpp"
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
    auto solve_(SystemT& sys, const VectorT& b) {

      double tolerance = this->getTolerance();
      double nonLinearTolerance = this->getNonLinearTolerance();

      core::defaultTypes::uint maxIterations = this->getMaxIterations();
      core::defaultTypes::uint maxNonLinearIterations =
	this->getMaxNonLinearIterations();

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
