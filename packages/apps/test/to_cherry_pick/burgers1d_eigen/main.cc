
//#include "ode_observer.hpp"
// app class
#include "apps_burgers1d_eigen.hpp"
// integrator
#include "integrators/ode_integrate_n_steps.hpp"
// steppers
#include "explicit_steppers/ode_explicit_euler_stepper.hpp"
//#include "step_methods/ode_explicit_runge_kutta4_stepper.hpp"


struct mysizer{
  using state_t = Eigen::VectorXd;
 static size_t getSize(const state_t & obj){
   return obj.size();
 };
 static void resize(state_t & obj, size_t newSize){
   obj.resize(newSize);
 };
};


int main(int argc, char *argv[])
{
  using state_t = apps::burgers1dEigen::state_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;

  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, 20);
  appObj.setup();

  // integrate in time starting from y0
  scalar_t final_time = 35.;
  scalar_t dt = 0.001;
  
  {
    //********************************************
    // EXPLICIT EULER
    //********************************************
    state_t y0 = appObj.getInitialState();
    //    snapshot_collector collectorObj;

    using stepper_t = ode::explicitEulerStepper<
      state_t, state_t, scalar_t, 
      model_eval_t, scalar_t, mysizer
      /*default = standard policy */>;
    stepper_t stepperObj(appObj);
    
    ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_time/dt );
    // std::cout << collectorObj.getCount() << std::endl;
    // //    collectorObj.print();

    std::cout << "Final solution " << std::endl;
    for (int i=0; i<y0.size(); ++i)
      std::cout << std::setprecision(14) << y0(i)  << " ";
    std::cout << std::endl;
  }
  // {
  //   //********************************************
  //   // EXPLICIT RUNGE KUTTA 4
  //   //********************************************
  //   state_t y0 = appObj.getInitialState();
  //   snapshot_collector collectorObj;

  //   using stepper_t = ode::explicitRungeKutta4Stepper<
  //     state_t, state_t, scalar_t, eigenStateResizer,
  //     model_eval_t, scalar_t /*default = standard policy */>;
  //   stepper_t stepperObj(appObj);

  //   ode::integrateNSteps( stepperObj, y0, collectorObj, 0.0, dt, final_time/dt );
  //   std::cout << collectorObj.getCount() << std::endl;
  //   //    collectorObj.print();

  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0(i)  << " ";
  //   std::cout << std::endl;
  // }
  
  return 0;
}
