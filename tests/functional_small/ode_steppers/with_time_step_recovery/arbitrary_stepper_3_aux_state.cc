
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp
{
  std::string & checkStr_;
  using time_type   = double;
  using state_type    = Eigen::VectorXd;
  using discrete_time_residual_type = state_type;
  using discrete_time_jacobian_type = Eigen::SparseMatrix<double>;

public:
  MyApp(std::string & strin) : checkStr_(strin){}

  state_type createState() const{
    state_type ss(3);
    ss.setConstant(0);
    return ss;
  }

  discrete_time_residual_type createDiscreteTimeResidual() const{
    discrete_time_residual_type R(3);
    return R;
  }

  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
    discrete_time_jacobian_type J(3,3);
    return J;
  }

  template <typename step_t>
  void discreteTimeJacobian(const step_t & step,
                            const time_type & evaltime,
                            const time_type & dt,
                            discrete_time_jacobian_type & J,
        const state_type & yn,
        const state_type & ynm1,
        const state_type & ynm2,
        const state_type & ynm3) const
  {
    // dummy, not used for this test
  }

  template <typename step_t, typename state_type>
  void discreteTimeResidual(const step_t & step,
				const time_type & evaltime,
				const time_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1,
				const state_type & ynm2,
				const state_type & ynm3) const
  {
    std::cout << "-----------------\n";
    std::cout << step << "\n";

    std::cout << "yn: " << " "
	      << yn[0] << " "
	      << yn[1] << " "
	      << yn[2] << "\n";
    std::cout << "y_n-1: " << " "
	      << ynm1[0] << " "
	      << ynm1[1] << " "
	      << ynm1[2] << "\n";
    std::cout << "y_n-2: " << " "
	      << ynm2[0] << " "
	      << ynm2[1] << " "
	      << ynm2[2] << "\n";
    std::cout << "y_n-3: " << " "
	      << ynm3[0] << " "
	      << ynm3[1] << " "
	      << ynm3[2] << "\n";

    R[0] = ynm1[0]+dt;
    R[1] = ynm1[1]+dt;
    R[2] = ynm1[2]+dt;

    if (step==1){
      // for ode step1, yn changes but the others done
      if( std::abs(ynm1[0]-1.) > 1e-13 or
	  std::abs(ynm1[1]-1.) > 1e-13 or
	  std::abs(ynm1[2]-1.) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm2[0]-0) > 1e-13 or
	  std::abs(ynm2[1]-0) > 1e-13 or
	  std::abs(ynm2[2]-0) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm3[0]-0) > 1e-13 or
	  std::abs(ynm3[1]-0) > 1e-13 or
	  std::abs(ynm3[2]-0) > 1e-13){
	checkStr_ = "FAILED";
      }
    }

    if (step==2){
      if( std::abs(ynm1[0]-1.2) > 1e-13 or
	  std::abs(ynm1[1]-1.4) > 1e-13 or
	  std::abs(ynm1[2]-1.6) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm2[0]-1) > 1e-13 or
	  std::abs(ynm2[1]-1) > 1e-13 or
	  std::abs(ynm2[2]-1) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm3[0]-0) > 1e-13 or
	  std::abs(ynm3[1]-0) > 1e-13 or
	  std::abs(ynm3[2]-0) > 1e-13){
	checkStr_ = "FAILED";
      }
    }

    if (step==2 and (dt==0.1)){
      throw pressio::eh::DiscreteTimeResidualFailureUnrecoverable();
    }

    if (step==3){
      if( std::abs(ynm1[0]-1.4) > 1e-13 or
	  std::abs(ynm1[1]-1.8) > 1e-13 or
	  std::abs(ynm1[2]-2.2) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm2[0]-1.2) > 1e-13 or
	  std::abs(ynm2[1]-1.4) > 1e-13 or
	  std::abs(ynm2[2]-1.6) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(ynm3[0]-1) > 1e-13 or
	  std::abs(ynm3[1]-1) > 1e-13 or
	  std::abs(ynm3[2]-1) > 1e-13){
	checkStr_ = "FAILED";
      }
    }

    if (step==3 and (dt==0.1 or dt==0.05)){
      throw pressio::eh::DiscreteTimeResidualFailureUnrecoverable();
    }

    if (step==3){
      // the only one to succeed at step=3 should be dt=0.025 and such that
      if( std::abs(R[0]-1.425) > 1e-13 or
	  std::abs(R[1]-1.825) > 1e-13 or
	  std::abs(R[2]-2.225) > 1e-13){
	checkStr_ = "FAILED";
      }
    }
  }
};

struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    // this does not have any meaning, but it mimics the
    // steps happening inside the real pressio nonlin solver
    std::cout << "\n";

    state_t R(3);
    for (int i=0; i<2; ++i)
    {
      std::cout << "solver iter = " << i << "\n";
      try{
	sys.residual(state, R);
	state(0) += 0.1;
	state(1) += 0.2;
	state(2) += 0.3;
      }
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
	throw ::pressio::eh::NonlinearSolveFailure();
      }
    }
  }
};


int main(int argc, char *argv[])
{
  using app_t		= MyApp;
  using state_t = typename app_t::state_type;

  auto dtManager =
    [](const ::pressio::ode::step_count_type & step,
       const typename app_t::time_type & time,
       typename app_t::time_type & dt,
       typename app_t::time_type & minDt,
       typename app_t::time_type & dtRedFactor)
    {
      dt = 0.1;
      minDt = 0.01;
      dtRedFactor=2.;
    };

  auto collector =
    [](const ::pressio::ode::step_count_type & step,
       const typename app_t::time_type & time,
       const state_t & y){};

  std::string checkStr= "PASSED";
  app_t appObj(checkStr);
  state_t y(3);
  pressio::ops::fill(y, 1);
  MyFakeSolver solver;

  auto stepperObj = pressio::ode::create_arbitrary_stepper<4>(appObj);
  pressio::ode::advance_to_target_time_with_time_step_recovery_and_observe
    (stepperObj, y, 0., 0.4, dtManager, collector, solver);

  std::cout << checkStr << std::endl;
  return 0;
}
