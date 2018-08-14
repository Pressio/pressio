
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_BASE_HPP

#include <iostream>
#include <memory>
#include <type_traits>

#include "system_traits.hpp"
#include "meta/core_meta_static_checks.hpp"
#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"


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
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        !details::system_traits<SystemT>::is_system,
        int
      >::type* = nullptr
    >
    void solve(const SystemT& system, const VectorT& x0) {
      std::cerr << "Error: the first argument to method solve must be of system type" << std::endl;
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
      typename SystemT,
      typename VectorT,
      typename std::enable_if<
        details::system_traits<SystemT>::is_system,
        int
      >::type* = nullptr
    >
    auto solve(const SystemT& system, const VectorT& x0) {
      return 0;
    }


    uint getMaxIterations() {
      return maxIters_;
    }

    double getTolerance() {
      return tolerance_;
    }


    void setMaxIterations(uint maxIters) {
      maxIters_ = maxIters;
    }


    void setTolerance(double tolerance) {
      tolerance_ = abs(tolerance);
    }


//  protected:

    NonLinearSolverBase() : maxIters_(100), tolerance_(1.0e-5) {}

  
  private:

    uint maxIters_;
    double tolerance_;

};

} //end namespace solvers

#endif
