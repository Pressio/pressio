
#include "pressio_rom_lspg.hpp"
#include "custom_decoder.hpp"

constexpr double dt = 0.5;

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

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const
  {
    return dense_matrix_type(N_, B.cols());
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     scalar_type time,
  		     dense_matrix_type & A) const
  {
    Eigen::MatrixXd J(N_, N_);
    for (auto i=0; i<N_; ++i){
      for (auto j=0; j<N_; j+=2) J(i,j) = 1.;
      for (auto j=1; j<N_; j+=2) J(i,j) = 2.;
    }
    A = J*B;
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    f.setConstant(1.);
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
    // std::cout << "--- " << step << std::endl;
    // std::cout << romState << std::endl;
    // std::cout << std::endl;
    // std::cout << fomStateConstRef_ << std::endl;

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
  /* Here we verify correctness of view fom stata
     for LSPG with continuous-time API.

     Let g(x) be the decoder where g(x) = Ax where A[:,:]=2

     - we do two steps, from t_0 -> t_1 -> t_2

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

  using odetag = pressio::ode::implicitmethods::Euler;
  // using problem_t  = pressio::rom::lspg::composeDefaultProblem<
  //   odetag, fom_t, decoder_t, rom_state_t>::type;
  // problem_t problem(appObj, decoderObj, romState, refState);
  auto problem = pressio::rom::lspg::createDefaultProblemUnsteady<odetag>(
    appObj, decoderObj, romState, refState);

  using solver_t = MyFakeSolver<rom_state_t,typename decoder_t::jacobian_type>;
  solver_t solver(fomSize, romSize);

  Observer Obs(checkStr, problem.currentFomStateCRef());

  pressio::ode::advanceNSteps(problem.stepperRef(),
			      romState, 0.0, dt, 2,
			      Obs, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
