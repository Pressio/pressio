
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP

#include <type_traits>


#include "solvers_nonlinear_base.hpp"


namespace solvers {
	
struct NonlinearSolvers {

  template <typename SolverT>
  static auto createSolver() {
  	return NonLinearSolverBase();
  }

};

} // end namespace solvers

#endif