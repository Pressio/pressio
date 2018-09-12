
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP

#include <iostream>
#include <memory>
#include <type_traits>

#include "core_ConfigDefs.hpp"
#include "system_traits.hpp"
#include "solvers_norms_fwd.hpp"
#include "solvers_linear_traits.hpp"
#include "solvers_meta_static_checks.hpp"


namespace solvers {


struct NonlinearSolvers; // Fwd declaration

/**
 * @brief Base class for nonlinear solver implemented through CRTP
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a nonlinear solver.
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class NonLinearSolvers.
 */
template <typename Derived>
class NonLinearSolverBase {

  public:

    /**
     * Raise an assertion as the non linear system supplied as input is invalid
     *
     * @param  system non linear system to be solved
     * @param  x0 initial solution guess
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
     * Solve the non linear system.
     *
     * @param system non linear system to be solved.
     * @param x0 initial solution guess.
     * @return solution vector.
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
      return this->underlying().template solve_<SolverT, PrecT, NormT>(sys, x0);
    }


    /**
     * @brief Solve the non linear system.
     *
     * @param system is the non linear system to be solved.
     * @param x0 is the solution hint.
     * @return solution vector.
     *
     * @section DESCRIPTION
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
      return this->template solve<SolverT, linear::DefaultPreconditioner, L2Norm, SystemT, VectorT>(sys, x0);
    }


    /**
     * Get the maximum number of iterations of the nonlinear solver.
     */
    core::defaultTypes::uint getMaxNonLinearIterations() {
      return maxNonLinearIterations_;
    }


    /**
     * Get the tolerance of the nonlinear solver.
     */
    double getNonLinearTolerance() {
      return nonLinearTolerance_;
    }


    /**
     * Set the maximum number of iterations of the nonlinear solver.
     *
     * @param maxNonLinearIterations maximum number of iterations of the nonlinear solver.
     */
    void setMaxNonLinearIterations(core::defaultTypes::uint maxNonLinearIterations) {
      maxNonLinearIterations_ = maxNonLinearIterations;
    }


    /**
     * Set the tolerance of the nonlinear solver.
     *
     * @param nonLinearTolerance tolerance of the nonlinear solver.
     */
    void setNonLinearTolerance(double nonLinearTolerance) {
      nonLinearTolerance_ = abs(nonLinearTolerance);
    }


  protected:

    NonLinearSolverBase() : maxNonLinearIterations_(100), nonLinearTolerance_(1.0e-5) {}


  private:

    Derived& underlying() {
      return static_cast<Derived&>(*this);
    }


    Derived const& underlying() const {
      return static_cast<Derived const&>(*this);
    }


  private:

    core::defaultTypes::uint maxNonLinearIterations_;
    double nonLinearTolerance_;

};

} //end namespace solvers

#endif
