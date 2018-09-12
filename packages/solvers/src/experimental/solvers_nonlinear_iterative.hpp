
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_HPP

#include <type_traits>

#include "core_ConfigDefs.hpp"
#include "solvers_nonlinear_base.hpp"
#include "solvers_nonlinear_factory.hpp"


namespace solvers {


struct NonLinearSolvers; // Fwd declaration

/**
 * Implements a non linear solver bases on a linear iterative solver.
 */
template <
  typename PolicyT
>
class NonLinearIterativeSolver
  : public NonLinearSolverBase<
      NonLinearIterativeSolver<
        PolicyT
      >
    >
{

  private:

    friend NonLinearSolvers;
    typedef NonLinearSolverBase<NonLinearIterativeSolver<PolicyT>> base_type;


  public:


    /**
     * Implements the solve method for a non linear solver.
     */
    template <
      typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT
    >
    auto solve_(const SystemT& sys, const VectorT& b) {

      double nonLinearTolerance = this->getNonLinearTolerance();
      core::defaultTypes::uint maxNonLinearIterations = this->getMaxNonLinearIterations();

      return PolicyT::template solve<SolverT, PrecT, NormT>(sys, b, maxIterations_, maxNonLinearIterations, tolerance_, nonLinearTolerance);
    }


    /**
     * Get the maximum number of iterations of the underlying linear iterative solver.
     */
    core::defaultTypes::uint getMaxIterations() {
      return maxIterations_;
    }


    /**
     * Get the tolerance of the underlying linear iterative solver.
     */
    double getTolerance() {
      return tolerance_;
    }


    /**
     * Set the maximum number of iterations of the underlying linear iterative solver.
     *
     * @param maxIterations maximum number of iterations of the underlying linear iterative solver.
     */
    void setMaxIterations(core::defaultTypes::uint maxIterations) {
      maxIterations_ = maxIterations;
    }

    /**
     * Set the tolerance of the underlying linear iterative solver.
     *
     * @param tolerance tolerance of the underlying linear iterative solver.
     */
    void setTolerance(double tolerance) {
      tolerance_ = abs(tolerance);
    }


  protected:

  	NonLinearIterativeSolver() : base_type(), maxIterations_(100), tolerance_(1.0e-5) {}


  private:

    core::defaultTypes::uint maxIterations_;
    double tolerance_;
};

} // end namespace solvers

#endif
