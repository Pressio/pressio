
#include "pressio_solvers.hpp"
#include "./problems/epetra_expon_data_fit_n11.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.5
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f
 */

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  std::string sentinel = "PASSED";

  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN11;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(problem.numUn_);
  x(0) = 1.3;
  x(1) = 6.5e-1;
  x(2) = 6.5e-1;
  x(3) = 7.0e-1;
  x(4) = 6.0e-1;
  x(5) = 3.0;
  x(6) = 5.0;
  x(7) = 7.0;
  x(8) = 2.0;
  x(9) = 4.5;
  x(10) = 5.5;

  /* NOTE: this problem with GN only works with line search */

  using hessian_t = containers::DenseMatrix<Eigen::MatrixXd>;

  // linear solver type
  using solver_tag  = solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  // using gn_t = pressio::solvers::nonlinear::composeGaussNewton_t<
  //   problem_t, linear_solver_t>;
  // gn_t GNSolver(problem, x, linSolver);
  auto GNSolver = pressio::solvers::nonlinear::createGaussNewton(problem,x,linSolver);

  GNSolver.setUpdatingCriterion(pressio::solvers::nonlinear::update::armijo);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << " ";
  std::cout << std::endl;

  for(auto i=0; i<problem.numUn_; i++)
    if( std::abs(x(i) - problem.trueS[i]) > 1e-6 ) sentinel="FAILED";

  std::cout << sentinel << std::endl;

  MPI_Finalize();
}