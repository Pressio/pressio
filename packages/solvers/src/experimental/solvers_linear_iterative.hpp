
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "solvers_linear_base.hpp"


namespace rompp {
namespace solvers {

// Forward declarations
struct LinearSolvers;


/**
 * @brief Class that implements a linear iterative solver
 */
template<
  typename SolverT,
  typename MatrixT
>
class LinearIterativeSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      LinearIterativeSolver<
        SolverT,
        MatrixT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, LinearIterativeSolver<SolverT, MatrixT>> base_type;


  public:

    LinearIterativeSolver(LinearIterativeSolver&& other) :
      base_type(std::move(other)), maxIters_(other.maxIters_), tolerance_(other.tolerance_) {}


    template <typename T>
    auto _solve(const T& b)
     -> decltype(this->getSolver()->solve(b)) {
      auto solver = this->getSolver();
      solver->setMaxIterations(this->getMaxIterations());
      solver->setTolerance(this->getTolerance());
      return solver->solve(b);
    }


    inline core::default_types::uint getMaxIterations() {
      return maxIters_;
    }


    void setMaxIterations(core::default_types::uint maxIters) {
      maxIters_ = maxIters;
    }


    inline double getTolerance() {
      return tolerance_;
    }


    void setTolerance(double tolerance) {
      tolerance_ = tolerance;
    }


  protected:

    LinearIterativeSolver() :
      base_type(), maxIters_(100), tolerance_(1.0e-6) {};


    LinearIterativeSolver(std::shared_ptr<SolverT> solver) :
      base_type(solver), maxIters_(100), tolerance_(1.0e-6) {};


  private:

    core::default_types::uint maxIters_;
    double tolerance_;
};

} //end namespace solvers
} //end namespace rompp

#endif
