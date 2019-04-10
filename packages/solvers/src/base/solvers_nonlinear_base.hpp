
#ifndef SOLVERS_NONLINEAR_BASE_HPP
#define SOLVERS_NONLINEAR_BASE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../../core/src/vector/core_vector_meta.hpp"

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

  template <typename system_t>
  void solve(const system_t & sys, typename system_t::state_type & x){
    static_assert( ::rompp::solvers::meta::is_legitimate_system_for_nonlinear_solver<
      system_t>::value,
		   "The system obj type is not valid for non-linear solver");

    this->underlying().solveImpl(sys, x);
  }

protected:
  NonLinearSolverBase()	 = default;
  ~NonLinearSolverBase() = default;
  NonLinearSolverBase(const NonLinearSolverBase &) = delete;

private:
  Derived& underlying(){ return static_cast<Derived&>(*this);}
  Derived const& underlying() const {
    return static_cast<Derived const&>(*this);
  }

};

}}//end namespace rompp::solvers
#endif
