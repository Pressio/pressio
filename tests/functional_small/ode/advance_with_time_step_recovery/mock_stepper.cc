
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

template<typename ode_state_type>
struct MyFakeStepper
{
  template<typename solver_type>
  void doStep(ode_state_type & odeState,
	      const double & t,
	      const double & dt,
	      const pressio::ode::step_count_type & step,
	      solver_type & solver)
  {
    if (step==3 and dt==0.1)
      throw pressio::eh::TimeStepFailure();
    if (step==5 and (dt==0.1 or dt==0.05))
      throw pressio::eh::TimeStepFailure();

    for (int i=0; i<odeState.size(); i++){ 
      odeState(i) += dt;
    }
  }
};

struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state){}
};


int main(int argc, char *argv[])
{
  using ode_state_t = Eigen::VectorXd;

  ode_state_t y(3);
  y(0) = 1.0; y(1) = 2.0; y(2) = 3.0;

  MyFakeStepper<ode_state_t> stepper;
  MyFakeSolver solver;
  std::string checkStr= "PASSED";

  auto dtManager = [](const ::pressio::ode::step_count_type & step,
		      const double & time,
		      double & dt,
		      double & minDt,
		      double & dtRedFactor)
		{
		  dt = 0.1;
		  minDt = 0.01;
		  dtRedFactor=2.;
		};

  auto collector =
    [&checkStr](const ::pressio::ode::step_count_type & step,
		const double & time,
		const ode_state_t & y)
    {
      if (step==1){
	if( std::abs(y(0)-1.1) > 1e-10 or
	    std::abs(y(1)-2.1) > 1e-10 or
	    std::abs(y(2)-3.1) > 1e-10)
	  checkStr = "FAILED";

	if (std::abs(time-0.1) > 1e-10) checkStr="FAILED";
      }
      if (step==3){
	if( std::abs(y(0)-1.25) > 1e-10 or
	    std::abs(y(1)-2.25) > 1e-10 or
	    std::abs(y(2)-3.25) > 1e-10)
	  checkStr = "FAILED";

	if (std::abs(time-0.25) > 1e-10) checkStr="FAILED";
      }
      if (step==5){
	if( std::abs(y(0)-1.375) > 1e-10 or
	    std::abs(y(1)-2.375) > 1e-10 or
	    std::abs(y(2)-3.375) > 1e-10)
	  checkStr = "FAILED";

	if (std::abs(time-0.375) > 1e-10) checkStr="FAILED";
      }
      // std::cout << step << " "
      // 	       << y(0) << " "
      // 	       << y(1) << std::endl;
    };

  pressio::ode::advance_to_target_time_with_time_step_recovery
    (stepper, y, 0., 0.5, dtManager, collector, solver);

  if( std::abs(y(0)-1.575) > 1e-10 or
      std::abs(y(1)-2.575) > 1e-10 or
      std::abs(y(2)-3.575) > 1e-10)
    checkStr = "FAILED";

  std::cout << *y.data() << std::endl;
  std::cout << checkStr << std::endl;
  return 0;
}
