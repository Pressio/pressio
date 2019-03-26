
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "APPS_BURGERS1D"
#include "../fom_gold_states.hpp"

constexpr double eps = 1e-12;

template <typename state_t>
struct observer{
  using matrix_t = Eigen::MatrixXd;

  size_t state_size_ {};
  matrix_t A_;
  size_t count_ {};
  state_t y0_;
  state_t yIncr_;

  observer(int N, int state_size, const state_t & y0)
    : state_size_(state_size),
      A_(state_size, N+1), //+1 to store also init cond
      y0_(y0),
      yIncr_(state_size){}

  void operator()(size_t step,
  		  double t,
  		  const state_t & y){
    yIncr_ = y - y0_;
    this->storeInColumn(yIncr_, count_);
    count_++;
  }

  void storeInColumn(const state_t & y, int j){
    for (auto i=0; i<y.size(); i++)
      A_(i,j) = y(i);
  }
  size_t getCount() const{ return count_;}
  void printAll() const{ std::cout << A_ << std::endl; }
};


template <typename T>
void checkSol(const T & y, const std::vector<double> & trueS){
  for (size_t i=0; i< trueS.size(); i++)
    assert(std::abs(y[i] - trueS[i]) < eps);
}

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::Burgers1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::residual_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  //-------------------------------
  // create app object
  constexpr int Ncell = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  app_t appObj(mu, Ncell);
  appObj.setup();
  auto & y0n = appObj.getInitialState();

  // types for ode
  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_rhs_t>;
  using ode_jac_t   = rompp::core::Matrix<app_jacob_t>;

  ode_state_t y(y0n);
  using stepper_t = rompp::ode::ImplicitStepper<
    rompp::ode::ImplicitEnum::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_solver_t = rompp::solvers::EigenIterative<
    rompp::solvers::linear::Bicgstab, ode_jac_t>;
  rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
  solverO.setTolerance(1e-13);
  solverO.setMaxIterations(200);

  // integrate in time
  scalar_t fint = 0.10;
  scalar_t dt = 0.01;
  auto Nsteps = static_cast<unsigned int>(fint/dt);

  // define observer
  observer<ode_state_t> Obs(Nsteps, Ncell, y);

  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, Obs, solverO);
  Obs.printAll();

  checkSol(y, rompp::apps::test::Burg1DtrueImpEulerN20t010);
  std::cout << std::setprecision(14) << *y.data();

  return 0;
}
