
#ifndef SOLVERS_ITERATIVE_BASE_HPP
#define SOLVERS_ITERATIVE_BASE_HPP

#include "../solvers_ConfigDefs.hpp"

namespace rompp { namespace solvers {

template<typename scalar_t>
struct IterativeBase
{
  using iteration_t = core::default_types::uint;

  /** Get the maximum number of iterations. */
  inline iteration_t getMaxIterations() {
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
  inline scalar_t getTolerance() {
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
  IterativeBase() = default;
  IterativeBase(const IterativeBase &) = delete;
  ~IterativeBase() = default;

protected:
  iteration_t maxIters_ = static_cast<iteration_t>(100);
  scalar_t tolerance_	= static_cast<scalar_t>(0.000001);
};


}} //end namespace rompp::solvers

#endif
