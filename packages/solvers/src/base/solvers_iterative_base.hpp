
#ifndef SOLVERS_ITERATIVE_BASE_HPP
#define SOLVERS_ITERATIVE_BASE_HPP

#include "../solvers_ConfigDefs.hpp"

namespace pressio { namespace solvers {

template<typename derived_t, typename scalar_t>
struct IterativeBase
{
  IterativeBase() = default;
  IterativeBase(const IterativeBase &) = delete;
  ~IterativeBase() = default;

  using iteration_t = containers::default_types::uint;

  /** Get the number of iterations performed. */
  iteration_t getNumIterationsExecuted() const {
    return static_cast<const derived_t &>(*this).getNumIterationsExecutedImpl();
  }

  /** Get the error after last step is executed */
  scalar_t getFinalError() const {
    return static_cast<const derived_t &>(*this).getFinalErrorImpl();
  }

  /** Get the maximum number of iterations. */
  iteration_t getMaxIterations() const {
    return maxIters_;
  }

  /**
   * Set the maximum number of iterations
   *
   * @param maxIters maximum number of iterations.
   */
  void setMaxIterations(iteration_t maxIters) {
    maxIters_ = maxIters;
  }

  /** Get the tolerance. */
  scalar_t getTolerance() const {
    return tolerance_;
  }

  /**
   * Set the tolerance of the solver.
   *
   * @param tolerance tolerance of the solver.
   */
  void setTolerance(scalar_t tolerance) {
    tolerance_ = tolerance;
  }

protected:
  iteration_t maxIters_ = static_cast<iteration_t>(100);
  scalar_t tolerance_	= static_cast<scalar_t>(0.000001);
};


}} //end namespace pressio::solvers

#endif
