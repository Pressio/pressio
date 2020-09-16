
#include "pressio_solvers.hpp"
#include "./problems/eigen_expon_data_fit_n2.hpp"


template <typename problem_t, typename state_t, typename solver>
void testC1(std::string & sentinel, 
            problem_t & problem, 
            state_t & x, 
            solver & GNSolver)
{
  GNSolver.setTolerance(1e-8);
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
  auto criterion = pressio::solvers::nonlinear::stop::whenGradientAbsoluteNormBelowTolerance;
  GNSolver.setStoppingCriterion(criterion);
  // 1e3 is chosen to test the convergence condition
  GNSolver.setTolerance(1e3);
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
  auto criterion = pressio::solvers::nonlinear::stop::whenGradientRelativeNormBelowTolerance;
  GNSolver.setStoppingCriterion(criterion);
  GNSolver.setTolerance(1e-5);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.41728158794844) > 1e-11 )or 
       (std::abs(x(1) - 0.26465375115797) > 1e-11) ){
    sentinel = "FAILED";
  }
}

int main()
{
  std::string sentinel = "PASSED";

  using namespace pressio;

  using problem_t   = solvers::test::EigenExpDataFitN2;
  using state_t	    = typename problem_t::state_type;
  using mat_t   = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // linear solver type
  using solver_tag	= solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate, linear_solver_t>;
  solver GNSolver(problem, x, linSolver);


  x[0] = 2.0; x[1] = 0.25;
  testC1(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x[0] = 2.0; x[1] = 0.25;
  testC2(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x[0] = 2.0; x[1] = 0.25;
  testC3(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  std::cout << std::setprecision(14) << *x.data() << std::endl;

}


