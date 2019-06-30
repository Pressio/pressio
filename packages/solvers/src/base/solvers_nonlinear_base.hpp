
#ifndef SOLVERS_NONLINEAR_BASE_HPP
#define SOLVERS_NONLINEAR_BASE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../../containers/src/vector/containers_vector_meta.hpp"

namespace rompp{ namespace solvers{

/**
 * @brief Base class for nonlinear solver
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a nonlinear solver.
 */
template <typename Derived>
struct NonLinearSolverBase {

  NonLinearSolverBase()	 = default;
  ~NonLinearSolverBase() = default;
  NonLinearSolverBase(const NonLinearSolverBase &) = delete;

  template <
    typename system_t,
    typename state_t
    >
  void solve(const system_t & sys, state_t & x){
    // static_assert( ::rompp::solvers::meta::is_legitimate_system_for_nonlinear_solver<
    //   system_t>::value,
    // 		   "The system obj type is not valid for non-linear solver");

    static_cast<Derived&>(*this).solveImpl(sys, x);
  }
};

}}//end namespace rompp::solvers
#endif
