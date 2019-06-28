
#ifdef HAVE_PYBIND11
#ifndef SOLVERS_PY_GAUSS_NEWTON_HPP
#define SOLVERS_PY_GAUSS_NEWTON_HPP

#include "../../solvers_forward_declarations.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_normal_eq_impl.hpp"

namespace rompp{ namespace solvers{ namespace iterative{

/*
 * for interfacing with python
 * for time being, no-line search
*/
template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename when_converged_t
  >
class PyGaussNewton<
  system_t, state_t, residual_t, jacobian_t,
  hessian_t, linear_solver_t, scalar_t, when_converged_t,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_array_pybind11<state_t>::value and
    ::rompp::algebra::meta::is_array_pybind11<residual_t>::value and
    ::rompp::algebra::meta::is_array_pybind11<jacobian_t>::value and
    ::rompp::algebra::meta::is_array_pybind11<hessian_t>::value
    >
  >
  : public NonLinearSolverBase<
  PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, when_converged_t
    >
  >,
    public IterativeBase<scalar_t>
{
  static_assert( mpl::not_same<system_t, pybind11::object>::value,
		 "The PyGaussNewton is supposed to have state/res/jac = pybind11::array_t, but the system_t should NOT be a pybind11::object to contain a Python object. If you need to use a GaussNewton solver on numpy strucrtures, just use one ready off the shelf. This PyGaussNewton is tailored so that you can use an Application class written in Python and use the time-integrators within Pressio as well as the ROMs methods inside Pressio.");

  using this_t = PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, when_converged_t>;

  // need to friend base (crpt) to grant it access to this private methods
  friend NonLinearSolverBase<this_t>;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<scalar_t>;

  using typename iterative_base_t::iteration_t;

  // this references a Python class where we search for all needed ops
  pybind11::object pythonOps_ = {};

  linear_solver_t linSolver_ = {};
  residual_t res_    = {};
  jacobian_t jac_    = {};
  hessian_t hess_    = {};
  state_t JTR_       = {};

  // delta is the correction
  state_t dy_     = {};

  // ytrail needed if/when line search is used
  state_t ytrial_    = {};

  // norms
  scalar_t normO_    = {};
  scalar_t normN_    = {};

  // dummy observer
  utils::impl::empty obsObj_ = {};

public:
  PyGaussNewton() = delete;
  PyGaussNewton(const PyGaussNewton &) = delete;
  ~PyGaussNewton() = default;

  PyGaussNewton(const system_t	 & system,
		const state_t	 & yState,
		linear_solver_t & linearSolverIn,
		pybind11::object ops)
    : pythonOps_{ops},
      linSolver_(linearSolverIn),
      res_(system.residual(yState)),
      jac_(system.jacobian(yState)),
      hess_(pythonOps_.attr("multiply1")(jac_, jac_, true)),
      JTR_{ state_t(const_cast<state_t &>(yState).request()) },
      dy_{ state_t(const_cast<state_t &>(yState).request()) },
      ytrial_{ state_t(const_cast<state_t &>(yState).request()) },
      normO_{0}, normN_{0},
      obsObj_{}
  {}

private:

  void solveImpl(const system_t & sys, state_t & y)
  {
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    // find out which norm to use
    using norm_t = typename impl::NormSelectorHelper<when_converged_t>::norm_t;

    // policy for evaluating the norm of a vector
    using norm_evaluator_t = impl::ComputeNormHelper<norm_t>;

    // policy to checking convergence
    using is_converged_t = impl::IsConvergedHelper<when_converged_t>;

    /* policy for computing line search factor (alpha) such that
     * the update is done with y = y + alpha dy
     * alpha = 1 default when user does not want line search
     */
    using line_search_tag  = ::rompp::solvers::iterative::gn::noLineSearch;
    using lsearch_helper_t = impl::LineSearchHelper<line_search_tag>;

    // alpha for taking steps
    scalar_t alpha = {};
    // storing residaul norm
    scalar_t normRes = {};
    scalar_t normRes0 = {};

    constexpr auto one = ::rompp::utils::constants::one<scalar_t>();
    constexpr auto negOne = ::rompp::utils::constants::negOne<scalar_t>();

#ifdef DEBUG_PRINT
    auto ss = std::cout.precision();
    std::cout.precision(14);
    auto reset = utils::io::reset();
    auto fmt1 = utils::io::cyan() + utils::io::underline();
    const auto convString = std::string(is_converged_t::description_);
    ::rompp::utils::io::print_stdout(fmt1, "PyGN normal eqns:", "criterion:",
				    convString, reset, "\n");
#endif

    // compute the initial norm of y (the state)
    norm_evaluator_t::evaluate(y, normO_);
    normN_ = {0};

    iteration_t iStep = 0;
    while (iStep++ <= iterative_base_t::maxIters_)
    {
#ifdef DEBUG_PRINT
      ::rompp::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::rompp::utils::io::print_stdout(fmt, "PyGN step", iStep,
				      utils::io::reset(), "\n");
#endif

      // residual norm for current state
      norm_evaluator_t::evaluate(res_, normRes);
      // std::cout << "residual" << std::endl;
      // pythonOps_.attr("myprint")(res_);

      // store initial residual norm
      if (iStep==1) normRes0 = normRes;

      // compute LHS: J^T*J
      pythonOps_.attr("multiply2")(jac_, jac_, hess_, true);
      // std::cout << "\n";
      // std::cout << "jacobian" << std::endl;
      // pythonOps_.attr("myprint")(jac_);
      // std::cout << "\n";
      // std::cout << "hessian" << std::endl;
      // pythonOps_.attr("myprint")(hess_);

#ifdef DEBUG_PRINT
      auto fmt1 = utils::io::magenta() + utils::io::bold();
      ::rompp::utils::io::print_stdout(fmt1, "GN_JSize =",
				      jac_.shape()[0], jac_.shape()[1],
				      "\n");
      ::rompp::utils::io::print_stdout(fmt1, "GN_HessianSize =",
				    hess_.shape()[0], hess_.shape()[1],
				    utils::io::reset(), "\n");
#endif

      // compute RHS: J^T*res
      pythonOps_.attr("multiply2")(jac_, res_, JTR_, true);
      pythonOps_.attr("scale")(JTR_, negOne);

      // solve normal equations
      linSolver_.attr("solve")(hess_, JTR_, dy_);

      // compute norm of the correction
      norm_evaluator_t::evaluate(dy_, normN_);

#ifdef DEBUG_PRINT
      ::rompp::utils::io::print_stdout(std::scientific,
				      "||R|| =", normRes,
				      "||R||(r) =", normRes/normRes0,
				      "||dy|| =", normN_,
				      utils::io::reset(),
				      "\n");
#endif

      // compute multiplicative factor if needed
      lsearch_helper_t::evaluate(alpha, y, ytrial_, dy_, res_, jac_, sys);

      // solution update: y = y + alpha*dy;
      ::rompp::algebra::ops::do_update(y, one, dy_, alpha);
      // std::cout << "PyGN PrintUpdateSol" << std::endl;
      // pythonOps_.attr("myprint")(dy_);
      // pythonOps_.attr("myprint")(y);
      // std::cout << std::endl;

      // check convergence (whatever method user decided)
      auto flag = is_converged_t::evaluate(y, dy_, normN_,
					   normRes, normRes0,
					   iStep,
					   iterative_base_t::maxIters_,
					   iterative_base_t::tolerance_);

      // if we have converged, exit
      if (flag) {
	break;
      }

      // store new norm into old variable
      normO_ = normN_;

      // compute residual and jacobian
      sys.residual(y, res_);
      sys.jacobian(y, jac_);

      }//loop

#if defined DEBUG_PRINT
    std::cout.precision(ss);
    ::rompp::utils::io::print_stdout(std::fixed);
#endif

  }//end solveImpl
};


}}}//end namespace rompp::solvers::iterative
#endif
#endif
