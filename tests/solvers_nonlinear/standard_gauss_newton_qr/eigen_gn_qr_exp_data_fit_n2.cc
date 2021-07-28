
#include "pressio_solvers.hpp"
#include "eigen_expon_data_fit_n2.hpp"

template <typename problem_t, typename state_t, typename solver>
void testC1(std::string & sentinel, 
            problem_t & problem, 
            state_t & x, 
            solver & GNSolver)
{
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.4173449278229) > 1e-9 )or 
       (std::abs(x(1) - 0.26464986197941) > 1e-9) ){
    sentinel = "FAILED";
  }  
}

template <typename problem_t, typename state_t, typename solver>
void testC2(std::string & sentinel, 
            problem_t & problem, 
            state_t & x, 
            solver & GNSolver)
{
  // setting 2 max iters so the solver will exit afte 2 steps
  GNSolver.setMaxIterations(2);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(16) << *x.data() << std::endl;
  if ( (std::abs(x(0) - 2.415361667771343) > 1e-9 )or 
       (std::abs(x(1) - 0.2648293802571118) > 1e-9) ){
    sentinel = "FAILED";
  }  
}


int main()
{
  std::string sentinel= "PASSED";
  
  using namespace pressio;

  using problem_t   = solvers::test::EigenExpDataFitN2<double>;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;
  state_w_t x(2); 
  x(0) = 2.0; x(1) = 0.25;

  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  auto GNSolver = pressio::nonlinearsolvers::createGaussNewtonQR(problem,x,qrSolver);

  x(0) = 2.0; x(1) = 0.25;
  testC1(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC2(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  // // print summary from timers
  // #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  // utils::TeuchosPerformanceMonitor::stackedTimersReportSerial();
  // #endif
  return 0;
}
