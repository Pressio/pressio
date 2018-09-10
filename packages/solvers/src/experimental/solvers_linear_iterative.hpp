
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP

#include "core_ConfigDefs.hpp"
#include "solvers_linear_base.hpp"


namespace solvers {

// Forward declarations
struct LinearSolvers;


/**
 * @brief Class that implements a linear iterative solver
 */
template<
  typename SolverT,
  typename MatrixT,
  typename PolicyT
>
class LinearIterativeSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      PolicyT,
      LinearIterativeSolver<
        SolverT,
        MatrixT,
        PolicyT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, PolicyT, LinearIterativeSolver<SolverT, MatrixT, PolicyT>> base_type;


  public:

    LinearIterativeSolver(LinearIterativeSolver&& other) :
      base_type(std::move(other)), maxIters_(other.maxIters_), tolerance_(other.tolerance_) {}


    template <typename T>
    auto _solve(const T& b) {
      auto solver = this->getSolver();
      return PolicyT::solve(solver, b, this->getMaxIterations(), this->getTolerance());
    }


    inline core::defaultTypes::uint getMaxIterations() {
      return maxIters_;
    }


    void setMaxIterations(core::defaultTypes::uint maxIters) {
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

    core::defaultTypes::uint maxIters_;
    double tolerance_;
};

} //end namespace solvers

#endif
