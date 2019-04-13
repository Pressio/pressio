
#ifndef SOLVERS_NEWTON_RAPHSON_HPP
#define SOLVERS_NEWTON_RAPHSON_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CORE_OPS"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"

namespace rompp{ namespace solvers{

template <typename scalar_t,
	  typename linear_solver_t,
	  ::rompp::mpl::enable_if_t<
	    ::rompp::mpl::is_default_constructible<linear_solver_t>::value
	    > * = nullptr
	  >
class NewtonRaphson
  : public NonLinearSolverBase<NewtonRaphson<scalar_t, linear_solver_t>>,
    public IterativeBase<scalar_t>{

  using this_t	= NewtonRaphson<scalar_t,linear_solver_t>;
  using base_t  = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normN_ = {};
  linear_solver_t linSolver_; // default construct this

public:
  NewtonRaphson() = default;
  NewtonRaphson(const NewtonRaphson &) = delete;
  ~NewtonRaphson() = default;

private:
  template <typename T>
  scalar_t normOfDifference(const T & v1, const T& v2) const{
    T dVec(v1 - v2);
    return ::rompp::core::ops::norm2(dVec);
  }

public:
  template <typename system_t>
  void solveImpl(const system_t & sys,
		 typename system_t::state_type & x){
#ifdef DEBUG_PRINT
    std::cout << " starting Newton-Raphson solve "
	      << " tol = " << this->tolerance_
	      << " maxIter = " << this->maxIters_
	      << std::endl;
#endif
    using state_t    = typename system_t::state_type;

    auto Residual = sys.residual(x);
    auto Jac = sys.jacobian(x);

    auto dx(x);
    state_t xOld = x;

    linSolver_.solve(Jac, Residual, dx);
    x -= dx;

    core::default_types::uint iStep = 0;
    while (iStep++ <= this->maxIters_ &&
           this->normOfDifference(xOld, x) > this->tolerance_)
    {
      xOld = x;
      sys.residual(x, Residual);
      sys.jacobian(x, Jac);

      linSolver_.resetLinearSystem(Jac);
      linSolver_.solve(Residual, dx);
      x -= dx;
    }
  }//solveImpl

};//class

}} //end namespace rompp::solvers
#endif


// #ifdef DEBUG_PRINT
//       std::cout << " GN step=" << iStep
//    << " norm(dx)= " << normN_
//    << std::endl;
// #endif

// #ifdef DEBUG_PRINT
//  std::cout << " GN converged! " <<
//      << " final norm(dx)= " << normN_
//      << std::endl;
// #endif
