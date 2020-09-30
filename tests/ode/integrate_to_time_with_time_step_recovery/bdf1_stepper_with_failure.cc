
#include "pressio_rom.hpp"

struct MyApp
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  mutable int counts = 0;

public:
  velocity_type createVelocity() const{
    velocity_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void velocity(const state_type & yn,
		const scalar_type& time,
		velocity_type & f) const
  {
    std::cout << "f: t=" << time << "\n";

    if (std::abs(time-0.1) < 1e-13)
    {
      ++counts;
      if (counts<=2)
	throw pressio::eh::velocity_failure_unrecoverable();
    }

    f[0] = 1.;
    f[1] = 1.;
    f[2] = 2.;
  }

  void jacobian(const state_type&, const scalar_type&, jacobian_type&) const{
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
    // this does not have any meaning, but it mimics the
    // steps happening inside the real pressio nonlin solver

    ++count_;

    state_t R(3);
    for (int i=0; i<2; ++i)
    {
      std::cout << i << "\n";
      try{
	std::cout << "s: state" << " "
		  << state[0] << " "
		  << state[1] << " "
		  << state[2] << std::endl;

	sys.residual(state, R);

	std::cout << "s: res" << " "
		  << R[0] << " "
		  << R[1] << " "
		  << R[2] << std::endl;

	state[0] += 0.1;
	state[1] += 0.2;
	state[2] += 0.3;
      }
      catch (::pressio::eh::residual_evaluation_failure_unrecoverable const &e){
	throw ::pressio::eh::nonlinear_solve_failure();
      }
    }

    // for count=4, if we get here it means that we have successfully
    // completed the second time step with failure, because it fails twice,
    // so the end of time step ==2 corresponds to count=4
    if (count_==4)
    {
      if( std::abs(state[0]-1.4) > 1e-13 or
    	  std::abs(state[1]-1.8) > 1e-13 or
    	  std::abs(state[2]-2.2) > 1e-13){
    	checkStr_ = "FAILED";
      }
      if( std::abs(R[0]-0.075) > 1e-13 or
    	  std::abs(R[1]-0.175) > 1e-13 or
    	  std::abs(R[2]-0.25) > 1e-13){
    	checkStr_ = "FAILED";
      }
    }

    if (count_==5)
    {
      if( std::abs(state[0]-1.6) > 1e-13 or
    	  std::abs(state[1]-2.2) > 1e-13 or
    	  std::abs(state[2]-2.8) > 1e-13){
    	checkStr_ = "FAILED";
      }
      if( std::abs(R[0]-0.0) > 1e-13 or
    	  std::abs(R[1]-0.1) > 1e-13 or
    	  std::abs(R[2]-0.1) > 1e-13){
    	checkStr_ = "FAILED";
      }
    }
  }
};

int main(int argc, char *argv[])
{
  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
  using nstate_t	= typename app_t::state_type;
  using nvelo_t	        = typename app_t::velocity_type;
  using njacobian_t	= typename app_t::jacobian_type;

  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nvelo_t>;
  using jac_t		= ::pressio::containers::SparseMatrix<njacobian_t>;

  auto dtManager =
    [](const ::pressio::ode::types::step_t & step,
       const sc_t & time,
       sc_t & dt,
       sc_t & minDt,
       sc_t & dtRedFactor)
    {
      dt = 0.1;
      minDt = 0.01;
      dtRedFactor=2.;
    };

  auto collector =
    [](const ::pressio::ode::types::step_t & step,
		const sc_t & time,
		const state_t & y)
    {};

  std::string checkStr= "PASSED";
  app_t appObj;
  MyFakeSolver solver(checkStr);

  state_t y(3);
  pressio::ops::fill(y, 1);

  using ode_tag = pressio::ode::implicitmethods::Euler;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  pressio::ode::advanceToTargetTimeWithTimeStepRecovery
    (stepperObj, y, 0., 0.4, solver, dtManager, collector);

  std::cout << checkStr << std::endl;
  return 0;
}
