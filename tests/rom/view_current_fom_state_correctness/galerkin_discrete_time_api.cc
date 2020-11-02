
#include "pressio_rom_galerkin.hpp"
#include "custom_decoder.hpp"

struct MyFakeApp
{
  const int N_ = {};

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using discrete_time_residual_type = residual_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N)
    : N_(N){}

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
  			    Args && ... states) const
  {
    discreteTimeResidualImpl(step, time,dt,R,
			     std::forward<Args>(states)...);
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
  				 Args && ... states) const
  {
    A = Eigen::MatrixXd::Identity(N_, N_) *B;
  }

  discrete_time_residual_type createDiscreteTimeResidual() const
  {
    discrete_time_residual_type R(N_);
    return R;
  }

  dense_matrix_type createApplyDiscreteTimeJacobianResult
  (const dense_matrix_type & B) const
  {
    dense_matrix_type A(N_, B.cols());
    return A;
  }

private:
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    R.setConstant(1);
  }
};

template<typename r_t, typename j_t>
struct MyFakeSolver
{
  int romSize_ = {};
  r_t R_;
  j_t J_;

  MyFakeSolver(int romSize)
    : romSize_(romSize),
      R_(romSize),
      J_(romSize, romSize)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    for (auto k=0; k<2; ++k)
    {
      // we call these just because we should but
      // here not needed because this is a fake test
      sys.residual(state, R_);
      sys.jacobian(state, J_);

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
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
    : sentinel_(sentinel), fomStateConstRef_(fomStateRef){}

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
    // so the ode state should be the new one,
    // and the fomState should be the one corresponding to
    // the romState at the last call to the solver
    if (step==1){
      trueRom.setConstant(3.);
      trueFom.setConstant(20.);
      if (!romState.isApprox(trueRom) or
    	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }

    // step=2 means that we just completed the step 1 from t_1->t_2
    // and the collector is called at the end of the step
    // so the ode state should be the new one,
    // and the fomState should be the one corresponding to
    // the romState at the last call to the solver
    if (step==2){
      trueRom.setConstant(5.);
      trueFom.setConstant(40.);
      if (!romState.isApprox(trueRom) or
    	  !fomStateConstRef_.isApprox(trueFom))
      	sentinel_ = "FAILED";
    }
  }
};

int main(int argc, char *argv[])
{
  /* Here we verify that given a Galerkin problem
     for the discrete-time API, the method to view the
     current fom state gives the correct result.

     * dt = 0.5, we do steps so we go from t_0->t_1->t_2
     * the fake solver, which increments solution vector by 1.
     * R=[1 1 ... 1]^T always for simplicity

     - g(x) = [ 2 2 2 2 ]  x
                2 2 2 2
		...
	      [ 2 2 2 2 ]

      The solver fakes such that it increments the state by 2
      every time the solver is called.
      Since we have two steps, the solver should be called twice.

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
  refState.setZero();

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 1.0);

  // using ode_tag = pressio::ode::implicitmethods::Arbitrary;
  using rom_jacobian_t = pressio::containers::DenseMatrix<Eigen::MatrixXd>;
  // using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  // using problem_t = pressio::rom::galerkin::composeDefaultProblem<
  //   ode_tag, fom_t, rom_state_t, rom_jacobian_t,
  //   decoder_t, stepper_order, stepper_n_states>::type;
  // problem_t Problem(appObj, decoderObj, romState, refState);
  auto Problem =
    pressio::rom::galerkin::createDefaultProblem<rom_jacobian_t,1,2>(appObj, decoderObj, romState, refState);

  Observer Obs(checkStr, Problem.currentFomStateCRef());

  auto & stepperObj = Problem.stepperRef();
  MyFakeSolver<rom_state_t, rom_jacobian_t> solver(romSize);
  pressio::ode::advanceNSteps(stepperObj, romState, 0.0, 0.5, 2, Obs, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
