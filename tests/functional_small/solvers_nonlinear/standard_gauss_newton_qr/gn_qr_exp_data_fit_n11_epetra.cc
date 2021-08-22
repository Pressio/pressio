
#include "pressio/solvers.hpp"
#include "epetra_expon_data_fit_n11.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.5
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f
 */

bool test1()
{
  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN11;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

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
  // GaussNewton solver
  using qr_solver_type = qr::QRSolver<mat_type, qr::TSQR>;
  qr_solver_type qrSolver;

  auto GNSolver = pressio::nonlinearsolvers::create_gauss_newtonQR(problem,x,qrSolver);

  GNSolver.setUpdatingCriterion(pressio::nonlinearsolvers::Update::Armijo);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<problem.numUn_; i++)
    if( std::abs(x(i) - problem.trueS[i]) > 1e-6 )
      return false;

  return true;
}


bool test2()
{
  using namespace pressio;

  using problem_t   = solvers::test::EpetraExpDataFitN11;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

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

  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;
  auto GNSolver = pressio::nonlinearsolvers::create_gauss_newtonQR(problem,x,qrSolver);

  GNSolver.setUpdatingCriterion(pressio::nonlinearsolvers::Update::Armijo);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<problem.numUn_; i++)
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
  if (!b1 or !b2) sentinel = "FAILED";

  std::cout << sentinel << std::endl;

  MPI_Finalize();
}
