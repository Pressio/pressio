#include "pressio_ode_implicit.hpp"
#include "pressio_apps.hpp"
std::string checkStr {"PASSED"};
template <typename state_t>
struct observer{
  std::ofstream myfile_;
  observer(): myfile_("solution.bin",  std::ios::out | std::ios::binary){}

  void operator()(size_t step,
  		  double t,
  		  const state_t & y){
    auto ydata = *y.data();
    for (int i=0;i<y.extent(0);i++){
      myfile_.write(reinterpret_cast<const char*>(&ydata(i)),sizeof(ydata(i)));
    }
    std::cout << "Time = " << t << std::endl; 
  }

  void closeFile(){
    myfile_.close();
  }
};

/*
   Regression test for the 2D Shallow Water Equations with Eigen data structures
   and Crank Nicolson time stepping. Solution will be written to a solution.bin 
   output file. 
*/

int main(int argc, char *argv[]){
  using scalar_t = double;
  using app_t	 = ::pressio::apps::swe2d;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;


  constexpr int nx = 8;
  constexpr int ny = 8;
  scalar_t params[3];
  params[0] = 9.8;
  params[1] = 0.125;
  params[2] = 0.25;
  constexpr scalar_t Lx = 5;
  constexpr scalar_t Ly = 5;
  app_t appObj(Lx,Ly,nx,ny,params);
  scalar_t t = 0;
  constexpr scalar_t et = 10.;
  constexpr scalar_t dt = 0.5;

  // types for ode
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;

  ode_state_t y(appObj.getGaussianIC(params[1]));
  using ode_tag = pressio::ode::implicitmethods::CrankNicolson;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;

  // Define observer
  observer<ode_state_t> Obs;

  auto NonLinSolver=
    pressio::solvers::nonlinear::createNewtonRaphson(stepperObj, y,linSolverObj);
  NonLinSolver.setTolerance(1e-11);
  // integrate in time
  auto Nsteps = static_cast<::pressio::ode::types::step_t>(et/dt);

  std::string filename = "solution.bin";
  std::ofstream myfile (filename,  std::ios::out | std::ios::binary);

  pressio::ode::advanceNSteps(stepperObj, y, t, dt, Nsteps,Obs, NonLinSolver);
  Obs.closeFile();
  constexpr double solNormGold = 8.1221307554237;
  auto solNorm = (*y.data()).norm();
  if (std::abs(solNorm - solNormGold) > 1e-12){
    checkStr = "Failed";
  }
  std::cout << checkStr << std::endl;
  return 0;
}
