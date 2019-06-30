#ifndef SOLVERS_NONLINEAR_ITERATIVE_HELPER_HPP
#define SOLVERS_NONLINEAR_ITERATIVE_HELPER_HPP

#include "../solvers_ConfigDefs.hpp"


namespace rompp {
namespace solvers {

class NonLinearIterativeSolverHelper {
  
  public:

    /**
     * Construct the helper object
     */
     NonLinearIterativeSolverHelper() : tolerance_(1.0e-5), maxIterations_(100) {}


    /**
     * Get the maximum number of iterations of the underlying linear iterative solver.
     */
    containers::default_types::uint getMaxIterations() {
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
    void setMaxIterations(containers::default_types::uint maxIterations) {
      maxIterations_ = maxIterations;
    }

    /**
     * Set the tolerance of the underlying linear iterative solver.
     *
     * @param tolerance tolerance of the underlying linear iterative solver.
     */
    void setTolerance(double tolerance) {
      tolerance_ = std::abs(tolerance);
    }


  private:

    double tolerance_;
    containers::default_types::uint maxIterations_;
};

} // end namespace solvers
} // end namespace rompp

#endif
