
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP

#include <iostream>
#include <memory>
#include <type_traits>

#include "core_ConfigDefs.hpp"
#include "system_traits.hpp"
#include "solvers_meta_static_checks.hpp"
// #include "matrix/core_matrix_traits_exp.hpp"
// #include "vector/core_vector_traits_exp.hpp"


namespace solvers {


/**
 * @brief Base class for nonlinear solver implemented through CRTP
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a nonlinear solver.
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class NonLinearSolvers.
 */
//template<
//  typename MatrixT,
//  typename PolicyT,
//  typename Derived
//>

struct NonlinearSolvers; // Fwd declaration

template <typename Derived>
class NonLinearSolverBase {

  public:

    /**
     * @brief  Raise an assertion as the non linear system supplied as input is invalid
     *
     * @param  system is the non linear system to be solved
     * @param  x0 is the solution vector
     */
    template <
      typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        !(
          details::system_traits<SystemT>::is_system &&
          solvers::meta::are_vector_matrix_compatible<
            VectorT,
            typename details::system_traits<SystemT>::matrix_type
          >::value
        ),
        int
      >::type* = nullptr
    >
    void solve(const SystemT& sys, const VectorT& x0) {
      std::cerr << "Error: either the nonlinear system or the solution hint is invalid." << std::endl;
      assert(0);
    }


    /**
     * @brief  Solve the non linear system
     *
     * @param  system is the non linear system to be solved
     * @param  x0 is the solution hint
     * @return Solution vector
     */
    template <
      typename SolverT,
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        details::system_traits<SystemT>::is_system &&
        solvers::meta::are_vector_matrix_compatible<
          VectorT,
          typename details::system_traits<SystemT>::matrix_type
        >::value,
        int
      >::type* = nullptr
    >
    auto solve(const SystemT& sys, const VectorT& x0) {
      return 0;
    }


    /**
     * @brief  Solve the non linear system
     *
     * @param  system is the non linear system to be solved
     * @param  x0 is the solution hint
     * @return Solution vector
     *
     * DESCRIPTION
     *
     * This version of solve takes a reduced set of meta-parameters
     * and forward the arguments tto the full solve method.
     */
    template <
      typename SolverT,
      typename SystemT,
      typename VectorT
    >
    auto solve(const SystemT& sys, const VectorT& x0) {
      return this->template solve<SolverT, void, void, SystemT, VectorT>(sys, x0);
    }


    core::defaultTypes::uint getMaxIterations() {
      return maxIters_;
    }

    double getTolerance() {
      return tolerance_;
    }


    void setMaxIterations(core::defaultTypes::uint maxIters) {
      maxIters_ = maxIters;
    }


    void setTolerance(double tolerance) {
      tolerance_ = abs(tolerance);
    }


  protected:

    NonLinearSolverBase() : maxIters_(100), tolerance_(1.0e-5) {}


  private:

    core::defaultTypes::uint maxIters_;
    double tolerance_;

};

} //end namespace solvers

#endif
