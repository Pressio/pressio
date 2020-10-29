
#include "pressio_rom_galerkin.hpp"
#include "custom_decoder.hpp"

struct MyFakeApp
{
  int N_ = {};

public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N)
    : N_(N){}

  velocity_type createVelocity() const{
    return state_type(N_);
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    f.setConstant(1.);
  }
};

struct Observer
{
  using fom_state_t = Eigen::VectorXd;
  using rom_state_t = Eigen::VectorXd;

  std::string & sentinel_;
  const fom_state_t & fomStateConstRef_;

  Observer(std::string & sentinel,
	   const fom_state_t & fomStateRef)
    : sentinel_(sentinel),
      fomStateConstRef_(fomStateRef){}

  void operator()(pressio::ode::types::step_t step,
		  double time,
		  const rom_state_t & romState)
  {
    std::cout << "--- " << step << std::endl;
    std::cout << romState << std::endl;
    std::cout << std::endl;
    std::cout << fomStateConstRef_ << std::endl;

    rom_state_t trueRom(romState.size());
    fom_state_t trueFom(fomStateConstRef_.size());

    // step=0 is because the collected is called BEFORE
    // starting the time loop
    if (step==0){
      trueRom.setConstant(1.);
      trueFom.setConstant(10.);
      if (!romState.isApprox(trueRom) or
	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }

    // step=1 means that we just completed the step 1 from t_0->t_1
    // and the collector is called at the end of the step
    // so the ode state should be the new one, but the
    // fomState should still be the one from before because
    // it has not been updated yet
    if (step==1){
      trueRom.setConstant(6.);
      trueFom.setConstant(10.);
      if (!romState.isApprox(trueRom) or
	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }

    // step=2 means that we just completed the step from t_1->t_2
    // and the collector is called at the end of the step
    // so the ode state should be the new one, but the
    // fomState should still be the one based on the previous state
    // before because it has not been updated yet
    if (step==2){
      trueRom.setConstant(11.);
      trueFom.setConstant(60.);
      if (!romState.isApprox(trueRom) or
	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }

    // step=3 means that we just completed the step from t_2->t_3
    // and the collector is called at the end of the step
    // so the ode state should be the new one, but the
    // fomState should still be the one based on the previous state
    // before because it has not been updated yet
    if (step==3){
      trueRom.setConstant(16.);
      trueFom.setConstant(110.);
      if (!romState.isApprox(trueRom) or
	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }
  }
};

int main(int argc, char *argv[])
{
  /* Here we verify that given a Galerkin problem
     for the continuous-time API, the method to view the
     current fom state gives the correct result.

	x_n+1 =  x_n + dt * phi^T f()

     To do this, we craft a test where:

     - f[:] = 1 always
     - x0[:] = 1
     - phi[:,:] = 1

     - we do 3 steps so from t_0->t_1->t_2->t_3
     and the mapping is g(x)=Ax where A[:,:]=2

     Based on this we know what fomState should
     be at every step.
     we use the observer to verify things behave correctly
   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 10;
  constexpr int romSize = 5;

  // app object
  fom_t appObj(fomSize);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);
  refState.setConstant(0);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 1.0);

  using ode_tag = pressio::ode::explicitmethods::Euler;
  // using problem_t  = pressio::rom::galerkin::composeDefaultProblem<
  //   ode_tag, fom_t, decoder_t, rom_state_t>::type;
  // problem_t galerkinProb(appObj, decoderObj, romState, refState);
  auto galerkinProb =
    pressio::rom::galerkin::createDefaultProblem<ode_tag>(appObj, decoderObj, romState, refState);

  Observer Obs(checkStr, galerkinProb.currentFomStateCRef());

  pressio::ode::advanceNSteps(galerkinProb.stepperRef(), romState,
   			      0.0, 0.5, 3, Obs);

  std::cout << checkStr <<  std::endl;
  return 0;
}
