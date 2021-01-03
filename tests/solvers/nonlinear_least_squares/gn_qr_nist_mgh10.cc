// Note: This problem is mgh10 from 
// https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
#include "pressio_solvers.hpp"
#include "./problems/eigen_nist_mgh10.hpp"


template <typename problem_t, typename state_t, typename solver>
void testC1(std::string & sentinel, 
            problem_t & problem, 
            state_t & x, 
            solver & GNSolver)
{
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 5.6096364710E-03) > 2.e-4) or 
       (std::abs(x(1) - 6.1813463463E+03) > 3.e+1) or 
       (std::abs(x(2) - 3.4522363462E+02) > 8.e-1) ){
    sentinel = "FAILED";
  }  
}

template <typename problem_t, typename state_t, typename solver>
void testC2(std::string & sentinel, 
            problem_t & problem, 
            state_t & x, 
            solver & GNSolver)
{
  // setting 2 max iters so the solver will exit after 20 steps
  GNSolver.setMaxIterations(20);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(16) << *x.data() << std::endl;
  if ( (std::abs(x(0) - 5.6096364710E-03) > 2.e-4) or 
       (std::abs(x(1) - 6.1813463463E+03) > 3.e+1) or 
       (std::abs(x(2) - 3.4522363462E+02) > 8.e-1) ){
    sentinel = "FAILED";
  }  
}


int main()
{
  std::string sentinel= "PASSED";
  
  using namespace pressio;

  using problem_t   = solvers::test::EigenNISTmgh10;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;
  state_w_t x(3); 
  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;

  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  auto GNSolver = pressio::solvers::nonlinear::createGaussNewtonQR(problem,x,qrSolver);

  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;
  testC1(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;
  testC2(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  // // print summary from timers
  // #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  // utils::TeuchosPerformanceMonitor::stackedTimersReportSerial();
  // #endif
}
