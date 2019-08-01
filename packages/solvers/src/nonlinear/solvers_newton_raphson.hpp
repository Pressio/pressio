
#ifndef SOLVERS_NEWTON_RAPHSON_HPP
#define SOLVERS_NEWTON_RAPHSON_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CONTAINERS_BASIC"
#include "../../../CONTAINERS_OPS"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#endif

namespace pressio{ namespace solvers{


// template <typename T, typename enable = void>
// struct LinearSolveDispatch;

// template <typename T>
// struct LinearSolveDispatch<
//   T,
//   mpl::enable_if_t<
//     mpl::not_same<T, pybind11::object>::value
//     >
//   >
// {
//   template <typename A_t, typename x_t, typename b_t>
//   static void solve(T & solver,
// 		    const A_t & A, x_t & x, const b_t & b){
//     solver.solve(A, b, x);
//   }
// };

// template <typename T>
// struct LinearSolveDispatch<
//   T,
//   mpl::enable_if_t<
//     mpl::is_same<T, pybind11::object>::value
//     >
//   >
// {
//   template <typename A_t, typename x_t, typename b_t>
//     static void solve(T & solver,
// 		      const A_t & A, x_t & x, const b_t & b){
//     //solver.solve(A, b, x);
//   }
// };
// //---------------------------------------------------------


template <
  typename scalar_t,
  typename linear_solver_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_default_constructible<linear_solver_t>::value
#ifdef HAVE_PYBIND11
    or std::is_same<linear_solver_t, pybind11::object>::value
#endif
    > * = nullptr
  >
class NewtonRaphson
  : public NonLinearSolverBase<NewtonRaphson<scalar_t, linear_solver_t>>,
    public IterativeBase<scalar_t>{

  using this_t	= NewtonRaphson<scalar_t,linear_solver_t>;
  using base_t  = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normN_ = {};
  linear_solver_t linSolver_;

public:
  NewtonRaphson() = default;
  NewtonRaphson(linear_solver_t & lsO) : linSolver_{lsO}{}
  NewtonRaphson(const NewtonRaphson &) = delete;
  ~NewtonRaphson() = default;

private:
  template <typename T>
  scalar_t normOfDifference(const T & v1, const T& v2) const{
    T dVec(v1 - v2);
    return ::pressio::containers::ops::norm2(dVec);
  }

public:
  template <
  typename system_t,
  typename state_t,
  mpl::enable_if_t<
    containers::meta::is_wrapper<state_t>::value
    > * = nullptr
  >
  void solveImpl(const system_t & sys, state_t & x)
  {
#ifdef DEBUG_PRINT
    std::cout << " starting Newton-Raphson solve "
	      << " tol = " << this->tolerance_
	      << " maxIter = " << this->maxIters_
	      << std::endl;
#endif

    auto dx(x);
    state_t xOld = x;
    normN_ = {0};
    containers::default_types::uint iStep = 0;
    std::cout.precision(15);
    while (++iStep <= this->maxIters_)
    {
#ifdef DEBUG_PRINT
      ::pressio::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::pressio::utils::io::print_stdout(fmt, "NewRaph step", iStep,
				      utils::io::reset(), "\n");
#endif

      xOld = x;
      auto Residual = sys.residual(x);
      auto Jac = sys.jacobian(x);

      linSolver_.solve(Jac, Residual, dx);
      x -= dx;
      normN_ =::pressio::containers::ops::norm2(dx);
      ::pressio::utils::io::print_stdout("norm(dx) =", normN_, "\n");
      if (normN_ < this->tolerance_)
      	break;
    }
  }//solveImpl


#ifdef HAVE_PYBIND11
  template <
    typename system_t,
    typename state_t,
    mpl::enable_if_t<
      containers::meta::is_cstyle_array_pybind11<state_t>::value and
      mpl::is_same<system_t, pybind11::object>::value and
      mpl::is_same<linear_solver_t, pybind11::object>::value
      > * = nullptr
    >
  void solveImpl(const system_t & sys, state_t & x)
  {
#ifdef DEBUG_PRINT
    std::cout << " starting Newton-Raphson solve "
	      << " tol = " << this->tolerance_
	      << " maxIter = " << this->maxIters_
	      << std::endl;
#endif

    auto dx = state_t(x.request());
    auto xOld = state_t(dx.request());
    auto Residual = sys.attr("residual1")(x);
    auto Jac = sys.attr("jacobian1")(x);
    normN_ = {0};
    containers::default_types::uint iStep = 0;
    std::cout.precision(15);
    while (++iStep <= this->maxIters_)
    {
#ifdef DEBUG_PRINT
      ::pressio::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::pressio::utils::io::print_stdout(fmt, "NewRaph step", iStep,
				      utils::io::reset(), "\n");
#endif

      xOld = x;
      sys.attr("residual2")(x, Residual);
      sys.attr("jacobian2")(x, Jac);

      linSolver_.attr("solve")(Jac, dx, Residual);
      for (auto i=0; i<x.size(); ++i)
	x.mutable_at(i) -= dx.at(i);

      normN_ =::pressio::containers::ops::norm2(dx);
      ::pressio::utils::io::print_stdout("norm(dx) =", normN_, "\n");
      if (normN_ < this->tolerance_)
      	break;
    }
  }//solveImpl
#endif


};//class

}} //end namespace pressio::solvers
#endif
