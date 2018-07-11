
#include "CORE_VECTOR"
#include "CORE_MATRIX"
#include "ODE_ALL"
#include "SOLVERS_EXP"
#include "apps_burgers1d_eigen.hpp"
#include "experimental/rom_galerkin_explicit_policy.hpp"
#include "observer.hpp"

struct mysizer{
 using state_t = core::vector<apps::burgers1dEigen::state_type>;
 static size_t getSize(state_t & obj){
   return obj.size();
 };
  static void matchSize(const state_t & src, state_t & obj){
    obj.resize(src.size());
 };
};

template<typename T>
void printSol(std::string mess, const T & y){
  std::cout << mess << std::endl;
  for (int i=0; i<y.size(); ++i)
    std::cout << std::setprecision(14) << y[i]  << " ";
  std::cout << std::endl;
}


int main(int argc, char *argv[])
{
  using native_state_t = apps::burgers1dEigen::state_type;
  using native_jac_t = apps::burgers1dEigen::jacobian_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;

  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, 10);
  appObj.setup();

  // wrap with core structures
  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using residual_t = state_t;
  native_state_t y0n = appObj.getInitialState();
  // init state vector
  state_t y(y0n);

  // policy for computing space residual using the target app
  using res_pol_t = rom::exp::romGalerkinResidualPolicy<
    state_t, residual_t, model_eval_t, scalar_t, mysizer>;
  res_pol_t resObj(y);

  //using stepper_t = ode::explicitEulerStepper<
  using stepper_t = ode::explicitRungeKutta4Stepper<
    state_t, residual_t, scalar_t, model_eval_t,
    scalar_t, mysizer, res_pol_t>;
  stepper_t stepperObj(appObj, resObj);

  // integration details
  scalar_t dt = 0.01;
  scalar_t final_t = dt*1;
  snapshot_collector collObj;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, final_t/dt, collObj);    
  printSol("", y);
       
  return 0;
}
