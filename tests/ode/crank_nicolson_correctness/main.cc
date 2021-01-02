
#include "pressio_rom.hpp"

struct MyApp
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  velocity_type createVelocity() const{ return velocity_type(3); }
  jacobian_type createJacobian() const{return jacobian_type(3,3);}

  void velocity(const state_type & y,
		const scalar_type& time,
		velocity_type & f) const
  {
    std::cout << "velo: t=" << time << "\n";
    f[0] = y(0)+time;
    f[1] = y(1)+time;
    f[2] = y(2)+time;
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
    ++count_;
    std::cout << "SOLVE count = "  << count_ << std::endl;

    state_t R(3);
    for (int i=0; i<2; ++i)
    {
      std::cout << "i = "  << i << std::endl;
      sys.residual(state, R);
      std::cout << "state = "  << *state.data() << std::endl;
      std::cout << "R = " << *R.data() << std::endl;

      if (count_==1 && i==0)
      {
	if( std::abs(state(0)-1.) > 1e-13 or
	    std::abs(state(1)-2.) > 1e-13 or
	    std::abs(state(2)-3.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(1.-1.-0.5*1.5*2.5-0.5*1.5*1.) ) > 1e-13 or
	    std::abs(R(1)-(2.-2.-0.5*1.5*3.5-0.5*1.5*2.) ) > 1e-13 or
	    std::abs(R(2)-(3.-3.-0.5*1.5*4.5-0.5*1.5*3.) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      if (count_==1 && i==1)
      {
	if( std::abs(state(0)-2.) > 1e-13 or
	    std::abs(state(1)-3.) > 1e-13 or
	    std::abs(state(2)-4.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(2.-1.-0.5*1.5*3.5-0.5*1.5*1.) ) > 1e-13 or
	    std::abs(R(1)-(3.-2.-0.5*1.5*4.5-0.5*1.5*2.) ) > 1e-13 or
	    std::abs(R(2)-(4.-3.-0.5*1.5*5.5-0.5*1.5*3.) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      if (count_==2 && i==0)
      {
	if( std::abs(state(0)-3.) > 1e-13 or
	    std::abs(state(1)-4.) > 1e-13 or
	    std::abs(state(2)-5.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(3.-3.-0.5*1.5*6. - 0.5*1.5*4.5) ) > 1e-13 or
	    std::abs(R(1)-(4.-4.-0.5*1.5*7. - 0.5*1.5*5.5) ) > 1e-13 or
	    std::abs(R(2)-(5.-5.-0.5*1.5*8. - 0.5*1.5*6.5) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      if (count_==2 && i==1)
      {
	if( std::abs(state(0)-4.) > 1e-13 or
	    std::abs(state(1)-5.) > 1e-13 or
	    std::abs(state(2)-6.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(4.-3.-0.5*1.5*7. - 0.5*1.5*4.5) ) > 1e-13 or
	    std::abs(R(1)-(5.-4.-0.5*1.5*8. - 0.5*1.5*5.5) ) > 1e-13 or
	    std::abs(R(2)-(6.-5.-0.5*1.5*9. - 0.5*1.5*6.5) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      if (count_==3 && i==0)
      {
	if( std::abs(state(0)-5.) > 1e-13 or
	    std::abs(state(1)-6.) > 1e-13 or
	    std::abs(state(2)-7.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(5.-5.-0.5*1.5*9.5  - 0.5*1.5*8.) ) > 1e-13 or
	    std::abs(R(1)-(6.-6.-0.5*1.5*10.5 - 0.5*1.5*9.) ) > 1e-13 or
	    std::abs(R(2)-(7.-7.-0.5*1.5*11.5 - 0.5*1.5*10.) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      if (count_==3 && i==1)
      {
	if( std::abs(state(0)-6.) > 1e-13 or
	    std::abs(state(1)-7.) > 1e-13 or
	    std::abs(state(2)-8.) > 1e-13){
	  checkStr_ = "FAILED";
	}
	if( std::abs(R(0)-(6.-5.-0.5*1.5*10.5 - 0.5*1.5*8.) ) > 1e-13 or
	    std::abs(R(1)-(7.-6.-0.5*1.5*11.5 - 0.5*1.5*9.) ) > 1e-13 or
	    std::abs(R(2)-(8.-7.-0.5*1.5*12.5 - 0.5*1.5*10.) ) > 1e-13)
	{
	  checkStr_ = "FAILED";
	}
      }

      state(0) += 1.;
      state(1) += 1;
      state(2) += 1;
    }
  }
};

int main(int argc, char *argv[])
{
  /*
    check numerics of crank-nicolson is doing the right thing numerically
    by manufacturing an integration problem where we use
    artificial numbers and we check insdie the "fake" solver that
    the operators are correct.

    we integrate: dy/dt = f(y,t)

    - start from t_0 = 0.0
    - y(t=0) = [1 2 3]
    - dt = 1.5
    - f=y+t where t is the time the f is evaluated at
    - the nonlinear solver fakes updating the solution by simply adding 1
    - at each time step we call the solver twice

    ********
    *step 1* from t_0 -> t_1
    ********
    - t_0 = 0
    - t_1 = 1.5
    this means we are at step = 0 and predicting solution at t_1
    so n = 0 and n+1 = 1

    - the first call to the solver should have:
    - R = y_1 - y_0 - 0.5*dt*f_1 - 0.5*dt*f_0, where f_0=f(y_0, t_0) and f_1=f(y_1,t_1)
    - where y_1 = y_0 because y_0 is the tentative guess
    so
    f_0 = y_0+0
    f_1 = y_1+t_1 = y_1+1.5

    - the second call to the solver should have:
    - R = y_1 - y_0 - 0.5*dt*f_1 - 0.5*dt*f_0, where f_0=f(y_0, t_0) and f_1=f(y_1,t_1)
    - y_1 is now the state coming out from the first solve attempt, so now y_1 = y_0+1
    so
    f_0 = y_0+0
    f_1 = y_1+t_1 = y_1+1.5


    ********
    *step 2* from t_1 -> t_2
    ********
    - t_1 = 1.5
    - t_2 = 3.
    this means we are at step = n = 1 and predicting solution at t_2
    so n = 1 and n+1 = 2

    at t_1 we should have y_1 = [3 4 5], this is the "result" from previous step

    - the first call to the solver should have:
    - R = y_2 - y_1 - 0.5*dt*f_2 - 0.5*dt*f_1, where f_1=f(y_1, t_1) and f_2=f(y_2,t_2)
    - where y_2 = y_1 because y_1 is the tentative guess
    so
    f_1 = y_1+t_1 = y_1+1.5
    f_2 = y_2+t_2 = y_2+3.0

    - the second call to the solver should have:
    - R = y_2 - y_1 - 0.5*dt*f_2 - 0.5*dt*f_1, where f_1=f(y_1, t_1) and f_2=f(y_2,t_2)
    - y_2 is now the state coming out from the first solve attempt, so now y_2 = y_1+1
    so
    f_1 = y_1+t_1 = y_1+1.5
    f_2 = y_2+t_2 = y_2+3.0

    ********
    *step 3* from t_2 -> t_3
    ********
    - t_1 = 3.
    - t_2 = 4.5
    this means we are at step = n = 2 and predicting solution at t_3
    so n = 2 and n+1 = 3

    at t_2 we should have y_2 = [5 6 7], this is the "result" from previous step

    - the first call to the solver should have:
    - R = y_3 - y_2 - 0.5*dt*f_3 - 0.5*dt*f_2,
    - where f_2=f(y_2, t_2) and f_3=f(y_3,t_3)
    - where y_3 = y_2 because y_2 is the tentative guess
    so
    f_1 = y_2+t_2 = y_2+3.0
    f_2 = y_3+t_3 = y_3+4.5

    - the second call to the solver should have:
    - R = y_3 - y_2 - 0.5*dt*f_3 - 0.5*dt*f_2,
    - y_3 is now the state coming out from the first solve attempt, so now y_3 = y_2+1
    so
    f_1 = y_2+t_2 = y_2+3.0
    f_2 = y_3+t_3 = y_3+4.5
   */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
  using nstate_t	= typename app_t::state_type;
  using nvelo_t	        = typename app_t::velocity_type;
  using njacobian_t	= typename app_t::jacobian_type;
  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nvelo_t>;
  using jac_t		= ::pressio::containers::SparseMatrix<njacobian_t>;

  std::string checkStr= "PASSED";
  app_t appObj;
  MyFakeSolver solver(checkStr);

  state_t y(3);
  y(0)=1.; y(1)=2.; y(2)=3.;

  using ode_tag = pressio::ode::implicitmethods::CrankNicolson;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);
  pressio::ode::advanceNSteps(stepperObj, y, 0., 1.5, 3, solver);
  std::cout << checkStr << std::endl;
  return 0;
}
