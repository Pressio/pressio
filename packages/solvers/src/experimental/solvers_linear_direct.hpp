
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_DIRECT_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_DIRECT_HPP

#include "solvers_linear_base.hpp"


namespace pressio{
namespace solvers{

// Forward declarations
struct LinearSolvers;


/**
 * Class that implements a linear direct solver.
 *
 * @tparam SolverT linear dense solver.
 * @tparam MatrixT matrix defining the linear system.
 * @tparam PolicyT policy that implements the solution algorithm.
 */
template<
  typename SolverT,
  typename MatrixT
//  typename PolicyT
>
class LinearDirectSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      LinearDirectSolver<
        SolverT,
        MatrixT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, LinearDirectSolver<SolverT, MatrixT>> base_type;


  public:

    LinearDirectSolver(LinearDirectSolver&& other) : base_type(std::move(other)) {}


    template <typename T>
    auto _solve(const T& b) 
      -> decltype(this->getSolver()->solve(b)) {
      return this->getSolver()->solve(b);
    }


  public:

    LinearDirectSolver() : base_type() {};

    LinearDirectSolver(std::shared_ptr<SolverT> solver) : base_type(solver) {};

};

} //end namespace solvers

}//end namespace pressio
#endif
