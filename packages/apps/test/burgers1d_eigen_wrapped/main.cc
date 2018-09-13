
#include "CORE_ALL"
#include "ODE_ALL"
// app class
#include "apps_burgers1d_eigen.hpp"
//#include "ode_observer.hpp"
//#include "../apps_helper_ode.hpp"

const std::vector<double> trueExpEuler =
  { 5.0209814000128,   5.044067908724,\
    5.0694601439541,  5.0973757621592,\
    5.1280505161248,  5.1617393082963,\
    5.1987172243105,  5.2392805237326,\
    5.2837475435221,  5.3324594086071,\
    5.3857802812742,  5.4440964817745,\
    5.5078129073313,  5.5773432783592,\
    5.6530870659136,  5.7353794504736,\
    5.8243903774842,  5.9199350492773,\
    6.0211454752168,  6.1259551255163};

const std::vector<double> trueExpRK4 =
  {
   5.0209814000127, 5.0440679087239,
   5.0694601439537, 5.0973757621556,
   5.1280505160954, 5.1617393080991,
   5.1987172232043, 5.2392805183854,
   5.2837475207746, 5.3324593220435,
   5.385779982776,  5.4440955393165,
   5.507810159816,  5.5773358328381,
   5.6530682119049, 5.7353346645434,
   5.8242903335432, 5.9197246664937,
   6.0207292250728, 6.1251822423728};

constexpr double eps = 1e-12;

template <typename T>
void checkSol(const T & y,
	      const std::vector<double> & trueS){
    assert(std::abs(y[0] - trueS[0]) < eps);
    assert(std::abs(y[1] - trueS[1]) < eps);
    assert(std::abs(y[2] - trueS[2]) < eps);
    assert(std::abs(y[3] - trueS[3]) < eps);
    assert(std::abs(y[4] - trueS[4]) < eps);

    assert(std::abs(y[5] - trueS[5]) < eps);
    assert(std::abs(y[6] - trueS[6]) < eps);
    assert(std::abs(y[7] - trueS[7]) < eps);
    assert(std::abs(y[8] - trueS[8]) < eps);
    assert(std::abs(y[9] - trueS[9]) < eps);

    assert(std::abs(y[10] - trueS[10]) < eps);
    assert(std::abs(y[11] - trueS[11]) < eps);
    assert(std::abs(y[12] - trueS[12]) < eps);
    assert(std::abs(y[13] - trueS[13]) < eps);
    assert(std::abs(y[14] - trueS[14]) < eps);

    assert(std::abs(y[15] - trueS[15]) < eps);
    assert(std::abs(y[16] - trueS[16]) < eps);
    assert(std::abs(y[17] - trueS[17]) < eps);
    assert(std::abs(y[18] - trueS[18]) < eps);
    assert(std::abs(y[19] - trueS[19]) < eps);
  
}//end method


int main(int argc, char *argv[]){
  
  //-------------------------------
  // define native types
  using app_state_t = apps::Burgers1dEigen::state_type;
  using app_space_residual_t = apps::Burgers1dEigen::space_residual_type;
  using scalar_t = apps::Burgers1dEigen::scalar_type;
  using target_app_t = apps::Burgers1dEigen;
  
  //-------------------------------
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  target_app_t appObj(mu, 20);
  appObj.setup();
  auto & y0n = appObj.getInitialState();
  auto & r0n = appObj.getInitialResidual();

  //-------------------------------
  // types for ode
  using ode_state_t = core::Vector<app_state_t>;
  using ode_res_t = core::Vector<app_space_residual_t>;
  
  ///////////////////
  // forward Euler
  ///////////////////
  {
    ode_state_t y(y0n);
    ode_res_t r(r0n);
    
    using stepper_t = ode::ExplicitEulerStepper<
      ode_state_t, ode_res_t, target_app_t>;
    stepper_t stepperObj(appObj, y, r);
    
    // integrate in time 
    scalar_t fint = 35;
    scalar_t dt = 0.01;
    ode::integrateNSteps(stepperObj, y, 0.0, dt,
			 static_cast<unsigned int>(fint/dt));
    //y.data()->Print(std::cout << std::setprecision(14));
    checkSol(y, trueExpEuler);
    std::cout << std::setprecision(13) << *y.data();// << std::endl;
  }

  ///////////////////
  // runge kutta 4th
  ///////////////////
  {
    ode_state_t y(y0n);
    ode_res_t r(r0n);

    using stepper_t = ode::ExplicitRungeKutta4Stepper<
      ode_state_t, ode_res_t, target_app_t>;
    stepper_t stepperObj(appObj, y, r);
  
    // integrate in time 
    scalar_t fint = 35;
    scalar_t dt = 0.01;
    ode::integrateNSteps(stepperObj, y, 0.0, dt, static_cast<unsigned int>(fint/dt));
    checkSol(y, trueExpRK4);
  }
  
  return 0;
}
