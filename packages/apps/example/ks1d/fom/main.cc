
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "APPS_KS1D"

#include <iostream>
#include <fstream>

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
    yIncr_ = y; // - y0_; // Just want to plot state for now
    this->storeInColumn(yIncr_, count_);
    count_++;
  }

  void storeInColumn(const state_t & y, int j){
    for (auto i=0; i<y.size(); i++)
      A_(i,j) = y(i);
  }

  void printAll() const
  {
	  std::ofstream myfile;
	  myfile.open("primal.dat");
	  myfile << A_ << std::endl;
	  myfile.close();
  }

  void printFinal() const
  {
	  std::ofstream myfile;
	  myfile.open("primal_final.dat");

	  for (auto i=0; i<y0_.size(); i++)
	  {
		  myfile << A_(i,count_-1) << std::endl;
	  }
	  myfile.close();
  }

};


int main(int argc, char *argv[]){
  using app_t		= rompp::apps::KS1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::residual_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  //-------------------------------
  // create app object
  constexpr int Nnodes = 127;
  Eigen::Vector3d mu(0.0, 0.0, 0.0);
  app_t appObj(mu, Nnodes);
  appObj.setup();
  auto & y0n = appObj.getInitialState();

  // types for ode
  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_rhs_t>;
  using ode_jac_t   = rompp::core::Matrix<app_jacob_t>;

  ode_state_t y(y0n);

  // define auxiliary stepper
  using aux_stepper_t = rompp::ode::ImplicitStepper<
    rompp::ode::ImplicitEnum::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  // nonimal stepper
  constexpr auto ode_case = rompp::ode::ImplicitEnum::BDF2;
  using stepper_t = rompp::ode::ImplicitStepper<
    ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

  // define solver
  using lin_solver_t = rompp::solvers::iterative::EigenIterative<
    rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
  solverO.setTolerance(1e-13);
  solverO.setMaxIterations(200);

  // integrate in time
  scalar_t fint = 1000.0;
  scalar_t dt = 0.1;
  auto Nsteps = static_cast<unsigned int>(fint/dt);

  // define observer
  observer<ode_state_t> Obs(Nsteps, Nnodes, y);

  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, Obs, solverO);
  Obs.printAll();
  Obs.printFinal();

  return 0;
}
