
#include "pressio_rom_lspg.hpp"
#include "custom_decoder.hpp"

constexpr double dt = 0.5;

struct MyFakeApp
{
  const int N_ = {};

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using discrete_time_residual_type = residual_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N) : N_(N){}

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

  template <typename step_t>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
				 const state_type & yn,
				 const state_type & ynm1) const
  {
    Eigen::MatrixXd J(N_, N_);
    for (auto i=0; i<N_; ++i){
      for (auto j=0; j<N_; j+=2) J(i,j) = 1.;
      for (auto j=1; j<N_; j+=2) J(i,j) = 2.;
    }
    A = B - dt*(J*B);
  }

  template <typename step_t>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
				 const state_type & yn,
				 const state_type & ynm1,
				 const state_type & ynm2) const
  {
    throw std::runtime_error("Not supposed to be called");
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
    // mimic the case where the velocity is computed
    // here the velocity, f, is always one
    Eigen::VectorXd f(N_); f.setConstant(1.);

    R = yn - ynm1 - dt*f;
  }
};


template<typename r_t, typename j_t>
struct MyFakeSolver
{
  int fomSize_ = {};
  int romSize_ = {};
  r_t R_;
  j_t J_;

  MyFakeSolver(int fomSize, int romSize)
    : fomSize_(fomSize),
      romSize_(romSize),
      R_(fomSize),
      J_(fomSize, romSize)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    Eigen::MatrixXd trueJ(fomSize_, romSize_);

    for (auto k=0; k<2; ++k)
    {
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

  }
};

int main(int argc, char *argv[])
{
  /* Here we verify correctness of
     viewing the currentFomState

     Let g(x) be the decoder, and let Jg be its jacobian.

     - we do 2 steps so from t_0->t_1->t_2

     - g(x) = [ 2 2 2 ..]  x
                2 2 2
		...
	      [ 2 2 2 ..]

     - the solver fakes a solution by incrementing x by 2 at every call

     - we use the observer to verify the fomState is correct
   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 7;
  constexpr int romSize = 5;

  // app object
  fom_t appObj(fomSize);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);
  refState.setConstant(0.);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 1.0);

  // using ode_tag = pressio::ode::implicitmethods::Arbitrary;
  // using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  // using problem_t  = pressio::rom::lspg::composeDefaultProblem<
    // ode_tag, fom_t, decoder_t, rom_state_t, stepper_order, stepper_n_states>::type;
  // problem_t problem(appObj, decoderObj, romState, refState);
  // 1=order of steppe, 2=tot # of states needed
  auto problem = pressio::rom::lspg::createDefaultProblemUnsteady<1,2>(
    appObj, decoderObj, romState, refState);


  const auto & currFomState = problem.currentFomStateCRef();

  // here, the fom state should be [10 10 ...]
  Eigen::VectorXd trueFomState(fomSize);
  trueFomState.setConstant(10.);
  if (!trueFomState.isApprox(currFomState)) checkStr = "FAILED";

  Observer Obs(checkStr, currFomState);

  MyFakeSolver<rom_state_t, typename decoder_t::jacobian_type> solver(fomSize, romSize);
  pressio::rom::lspg::solveNSequentialMinimizations(problem, romState,
  			      0.0, dt, 2, Obs, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
