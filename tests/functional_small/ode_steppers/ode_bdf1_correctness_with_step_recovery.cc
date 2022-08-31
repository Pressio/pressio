
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{ return state_type(3); }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void rightHandSide(const state_type & /*unused*/,
		const independent_variable_type& evaltime,
		right_hand_side_type & f) const
  {
    std::cout << "f: t=" << evaltime << "\n";

    if (std::abs(evaltime-0.2) < 1e-13)
    {
      throw pressio::eh::VelocityFailureUnrecoverable();
    }

    f[0] = 1.;
    f[1] = 1.;
    f[2] = 2.;
  }

  void jacobian(const state_type&,
		const independent_variable_type&,
		jacobian_type&) const{
    // not used by the test
  }
};

struct MyFakeSolver
{
  std::string & checkStr_;
  int count_={};

  MyFakeSolver(std::string & strin) : checkStr_(strin){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++count_;
    std::cout << "SOLVE count = "  << count_ << std::endl;


    auto R = sys.createResidual();
    auto J = sys.createJacobian();
    for (int i=0; i<2; ++i)
    {
      std::cout << i << "\n";
      try{
	std::cout << "s: state" << " "
		  << state(0) << " "
		  << state(1) << " "
		  << state(2) << std::endl;

	sys.residualAndJacobian(state, R, J, true);

	std::cout << "s: res" << " "
		  << R(0) << " "
		  << R(1) << " "
		  << R(2) << std::endl;

	state(0) += 0.1;
	state(1) += 0.2;
	state(2) += 0.3;
      }
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
	throw ::pressio::eh::NonlinearSolveFailure();
      }
    }

    if (count_==1)
    {
      if( std::abs(state(0)-1.2) > 1e-13 or
	  std::abs(state(1)-1.4) > 1e-13 or
	  std::abs(state(2)-1.6) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(R(0)-0.0) > 1e-13 or
	  std::abs(R(1)-0.1) > 1e-13 or
	  std::abs(R(2)-0.1) > 1e-13){
	checkStr_ = "FAILED";
      }
    }

    // at count==2 we have exception thrown

    if (count_==3)
    {
      if( std::abs(state(0)-1.4) > 1e-13 or
	  std::abs(state(1)-1.8) > 1e-13 or
	  std::abs(state(2)-2.2) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(R(0)-(1.3-1.2-0.05*1.)) > 1e-13 or
	  std::abs(R(1)-(1.6-1.4-0.05*1.)) > 1e-13 or
	  std::abs(R(2)-(1.9-1.6-0.05*2.)) > 1e-13)
      {
	checkStr_ = "FAILED";
      }
    }

    if (count_==4)
    {
      if( std::abs(state(0)-1.6) > 1e-13 or
	  std::abs(state(1)-2.2) > 1e-13 or
	  std::abs(state(2)-2.8) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(R(0)-0.0) > 1e-13 or
	  std::abs(R(1)-0.1) > 1e-13 or
	  std::abs(R(2)-0.1) > 1e-13){
	checkStr_ = "FAILED";
      }
    }
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using app_t		= MyApp;
  using state_t	= typename app_t::state_type;

  auto dtManager = [](pressio::ode::StepCount /*unused*/,
		      pressio::ode::StepStartAt<double> /*unused*/,
		      pressio::ode::StepSize<double> & dt,
		      pressio::ode::StepSizeMinAllowedValue<double> & minDt,
		      pressio::ode::StepSizeScalingFactor<double> & dtRedFactor)
  {
    dt = 0.1;
    minDt = 0.01;
    dtRedFactor=2.;
  };

  std::string checkStr= "PASSED";
  app_t appObj;
  MyFakeSolver solver(checkStr);

  state_t y(3);
  pressio::ops::fill(y, 1);

  auto stepperObj = pressio::ode::create_bdf1_stepper(appObj);
  pressio::ode::advance_to_target_point_with_step_recovery
    (stepperObj, y, 0., 0.4, dtManager, solver);

  std::cout << checkStr << std::endl;

  pressio::log::finalize();
  return 0;
}
