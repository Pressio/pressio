
#ifndef SOLVERS_NEWTON_RAPHSON_HPP
#define SOLVERS_NEWTON_RAPHSON_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CORE_OPS"

namespace rompp{ namespace solvers{


template <typename scalar_t,
	  typename linear_solver_t,
	  core::meta::enable_if_t<
	    core::meta::is_default_constructible<linear_solver_t>::value
	    > * = nullptr
	  >
class NewtonRaphson {

  using uint_t = core::default_types::uint;

private:
  uint_t maxNonLinearIterations_     = 500;
  scalar_t nonLinearTolerance_       = 1e-6;
  scalar_t normN_ 		               = {};
  linear_solver_t linSolver_; // default construct this

public:
  NewtonRaphson() = default;
  NewtonRaphson(const NewtonRaphson &) = delete;
  ~NewtonRaphson() = default;

public:
  void setMaxNonLinearIterations(uint_t maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(scalar_t nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }

  template <typename T>
  scalar_t normOfDifference(const T & v1, const T& v2) const{
    T dVec(v1 - v2);
    return ::rompp::core::ops::norm2(dVec);
  }



  template <typename SystemT,
	    typename VectorT,
	    core::meta::enable_if_t<
	      core::meta::is_core_vector_wrapper<VectorT>::value and
	      std::is_same<VectorT, typename SystemT::state_type >::value
	      > * =nullptr
	    >
  void solve(const SystemT& sys, VectorT& x)
  {
#ifdef DEBUG_PRINT
    std::cout << " starting Newton-Raphson solve "
	      << " tol = " << nonLinearTolerance_
	      << " maxIter = " << maxNonLinearIterations_
	      << std::endl;
#endif

    auto Residual = sys.residual(x);
    auto Jac = sys.jacobian(x);

    auto dx(x);
    VectorT xOld = x;

    linSolver_.solve(Jac, Residual, dx);
    x -= dx;

    core::default_types::uint iStep = 1;
    while (iStep++ <= maxNonLinearIterations_ &&
           normOfDifference(xOld, x) > nonLinearTolerance_)
    {
      xOld = x;
      sys.residual(x, Residual);
      sys.jacobian(x, Jac);

      linSolver_.resetLinearSystem(Jac);
      linSolver_.solve(Residual, dx);
      x -= dx;
    }

  }//solve

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
