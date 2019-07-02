
#ifndef SOLVERS_IMPL_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP
#define SOLVERS_IMPL_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../meta/solvers_is_legitimate_residual_observer_each_solver_step.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template<typename observer_t,
	 typename residual_t,
	 typename = void>
struct ResidualObserverEachSolverStep{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid,
		       int step) {
    // no op
  }
};


template<typename observer_t,
	 typename residual_t>
struct ResidualObserverEachSolverStep<
  observer_t, residual_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::solvers::meta::is_legitimate_residual_observer_each_solver_step<
      observer_t, residual_t
     >::value
    >
  >{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid,
		       int gn_step) {
    obs->observeResidualEachGNStep(gn_step, resid);
  }
};


}}}} //end namespace pressio::solvers::iterative::impl
#endif
