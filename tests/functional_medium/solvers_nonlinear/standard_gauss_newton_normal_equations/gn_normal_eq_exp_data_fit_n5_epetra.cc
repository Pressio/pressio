
#include "pressio_solvers.hpp"
#include "epetra_expon_data_fit_n5.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.4
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f
 */

bool test1()
{
  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN5;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(problem.numUn_);
  x(0) = 0.5;  x(1) = 1.5;
  x(2) = -1.0; x(3) = 0.01;
  x(4) = 0.02;

  using hessian_t = Eigen::MatrixXd;

  // linear solver type
  using solver_tag  = linearsolvers::iterative::LSCG;
  using linear_solver_t = linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  auto GNSolver = pressio::nonlinearsolvers::create_gauss_newton(problem,x,linSolver);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    if( std::abs(x(i) - problem.trueS[i]) > 1e-6 )
      return false;

  return true;
}

bool test2()
{
  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN5;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(5);
  x(0) = 0.5;  x(1) = 1.5;
  x(2) = -1.0; x(3) = 0.01;
  x(4) = 0.02;

  using hessian_t = Eigen::MatrixXd;

  // linear solver type
  using solver_tag  = linearsolvers::iterative::LSCG;
  using linear_solver_t = linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  auto GNSolver = pressio::nonlinearsolvers::create_gauss_newton(problem,x,linSolver);

  GNSolver.setUpdatingCriterion(pressio::nonlinearsolvers::Update::armijo);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    if( std::abs(x(i) - problem.trueS[i]) > 1e-6 ) return false;

  return true;
}

bool test3()
{
  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using hessian_t   = Eigen::MatrixXd;

  problem_t problem;
  state_t x(5);
  x(0) = 0.5;  x(1) = 1.5;
  x(2) = -1.0; x(3) = 0.01;
  x(4) = 0.02;

  using hessian_t = Eigen::MatrixXd;

  // linear solver type
  using solver_tag  = linearsolvers::iterative::LSCG;
  using linear_solver_t = linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  auto GNSolver = pressio::nonlinearsolvers::create_gauss_newton(problem,x,linSolver);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    if( std::abs(x(i) - problem.trueS[i]) > 1e-6 )
      return false;

  return true;
}


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  std::string sentinel = "PASSED";
  auto b1 = test1();
  auto b2 = test2();
  auto b3 = test3();
  if (!b1 or !b2 or !b3) sentinel = "FAILED";

  std::cout << sentinel << std::endl;

  MPI_Finalize();
  return 0;
}