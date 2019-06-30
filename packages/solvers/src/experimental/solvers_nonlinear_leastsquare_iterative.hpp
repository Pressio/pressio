#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP

#include <type_traits>
#include "solvers_nonlinear_base.hpp"
#include "solvers_nonlinear_iterative_helper.hpp"
#include "../solvers_ConfigDefs.hpp"


namespace rompp {
namespace solvers {

class NonLinearSolvers; // Fwd declaration

template <
  typename PolicyT,
  typename LSolverT
>
class NonLinearLeastSquareIterativeSolver
  : public NonLinearIterativeSolverHelper,
    public NonLinearSolverBase<
      NonLinearLeastSquareIterativeSolver<
        PolicyT,
        LSolverT
      >
    >
{

  private:

    friend NonLinearSolvers;
    typedef NonLinearSolverBase<NonLinearLeastSquareIterativeSolver<PolicyT, LSolverT>> base_type;


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
        PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b, 1, 1, 1e-2, 1e-2, 1.)
      ){

      double tolerance = this->getTolerance();
      double nonLinearTolerance = this->getNonLinearTolerance();
      containers::default_types::uint maxIterations = this->getMaxIterations();
      containers::default_types::uint maxNonLinearIterations = this->getMaxNonLinearIterations();

      return PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b, maxIterations, 
                maxNonLinearIterations, tolerance, nonLinearTolerance, lambda_);
    }


    /**
     * Get the value of lambda used in nonlinear least-square algorithms.
     */
    double getLambda() {
      return lambda_;
    }


  protected:

    NonLinearLeastSquareIterativeSolver() : NonLinearIterativeSolverHelper(), base_type(), lambda_(1.0) {}


  private:

    double lambda_;
};

} // end namespace solvers
} // end namespace rompp

#endif
