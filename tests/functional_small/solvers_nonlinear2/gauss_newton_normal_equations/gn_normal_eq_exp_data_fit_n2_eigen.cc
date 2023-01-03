
#include "pressio/solvers.hpp"
#include "eigen_expon_data_fit_n2.hpp"

template <typename problem_t, typename state_t, typename solver>
void testC1(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  GNSolver.stopperHandle().setTolerance(1e-8);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.4173449278229) > 1e-7 )or
       (std::abs(x(1) - 0.26464986197941) > 1e-7) ){
    sentinel = "FAILED";
  }
}

template <typename problem_t, typename state_t, typename solver>
void testC2(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  auto criterion = pressio::nonlinearsolvers::Stop::WhenAbsolutel2NormOfGradientBelowTolerance;
  // 1e3 is chosen to test the convergence condition
  GNSolver.stopperHandle().setCriterion(criterion, 1e3);
  GNSolver.solve(problem, x);

  if ( (std::abs(x(0) - 2.4153884777105201) > 1e-11 )or
       (std::abs(x(1) - 0.26879930127342189) > 1e-11) ){
    sentinel = "FAILED";
  }
}

template <typename problem_t, typename state_t, typename solver>
void testC3(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  auto criterion = pressio::nonlinearsolvers::Stop::WhenRelativel2NormOfGradientBelowTolerance;
  GNSolver.stopperHandle().setCriterion(criterion, 1e-5);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.41728158794844) > 1e-11 )or
       (std::abs(x(1) - 0.26465375115797) > 1e-11) ){
    sentinel = "FAILED";
  }
}


int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string sentinel = "PASSED";
  using namespace pressio;
  using problem_t   = solvers::test::EigenExpDataFitN2<double>;
  using state_t	    = typename problem_t::state_type;
  using mat_t   = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x(0) = 2.0; x(1) = 0.25;

  // linear solver type
  using solver_tag	= linearsolvers::iterative::LSCG;
  using linear_solver_t = linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  auto GNSolver = nonlinearsolvers::create_gauss_newton(problem, linSolver);

  x(0) = 2.0; x(1) = 0.25;
  testC1(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC2(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC3(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  std::cout << std::setprecision(14) << x << std::endl;
  pressio::log::finalize();
  return 0;
}
