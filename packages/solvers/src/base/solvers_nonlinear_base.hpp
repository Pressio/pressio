
#ifndef SOLVERS_NONLINEAR_BASE_HPP
#define SOLVERS_NONLINEAR_BASE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"

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

  template <
    typename system_t,
    core::meta::enable_if_t<
      core::meta::is_core_vector_wrapper<
	typename system_t::state_type>::value
      > * =nullptr
    >
  void solve(const system_t & sys, typename system_t::state_type & x){
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
