
#include "pressio_rom.hpp"
#include "custom_decoder.hpp"

struct MyFakeApp
{
  int N_ = {};

public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N) : N_(N){}

  residual_type createResidual() const{
    return state_type(N_);
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const
  {
    return dense_matrix_type(N_, B.cols());
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     dense_matrix_type & A) const
  {
    // Eigen::MatrixXd J(N_, N_);
    // for (auto i=0; i<N_; ++i){
    //   for (auto j=0; j<N_; j+=2) J(i,j) = 1.;
    //   for (auto j=1; j<N_; j+=2) J(i,j) = 2.;
    // }
    // A = J*B;
  }

  void residual(const state_type & state,
		residual_type & r) const
  {
    r.setConstant(1.);
  }
};

template<typename r_t, typename j_t>
struct MyFakeSolver
{
  using fom_state_t = Eigen::VectorXd;

  int fomSize_ = {};
  int romSize_ = {};
  int callCounter_ = 0;
  r_t R_;
  j_t J_;
  std::string & checkString_;
  const fom_state_t & fomStateConstRef_;

  MyFakeSolver(int fomSize,
	       int romSize,
	       std::string & checkString,
	       const fom_state_t & fomStateRef)
    : fomSize_(fomSize),
      romSize_(romSize),
      R_(fomSize),
      J_(fomSize, romSize),
      checkString_(checkString),
      fomStateConstRef_(fomStateRef)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++callCounter_;

    for (auto k=0; k<2; ++k)
    {
      std::cout << "Solver call " << callCounter_ << " " << k << std::endl;

      sys.residual(state, R_);
      sys.jacobian(state, J_);

      std::cout << *state.data() << std::endl;
      // std::cout << *R_.data() << std::endl;
      // std::cout << *J_.data() << std::endl;

      // at first call to solver, the fom state should
      // be [8 8 8 ...]
      if (callCounter_==1 and k==0)
      {
      	Eigen::VectorXd trueFomState(fomSize_);
      	trueFomState.setConstant(8.);
      	if (!trueFomState.isApprox(fomStateConstRef_)) checkString_ = "FAILED";
      }

      // at first call to solver, the fom state should
      // be [16 16 16 ...]
      if (callCounter_==1 and k==1)
      {
      	Eigen::VectorXd trueFomState(fomSize_);
      	trueFomState.setConstant(16.);
      	if (!trueFomState.isApprox(fomStateConstRef_)) checkString_ = "FAILED";
      }

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};

int main(int argc, char *argv[])
{
  /*
    Test to check fomState view behaves correctly.

    The mapping always computes g(x) = Ax, where A[:,:]=2
    the steady system always return f()=1

    the solver fakes a solution by adding 2 to the state

     we use the solver to verify that the fomState is correct
   */
  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 7;
  constexpr int romSize = 4;

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

  // using problem_type = typename pressio::rom::lspg::composeDefaultProblem<
  //     fom_t, decoder_t, rom_state_t>::type;
  // problem_type problem(appObj, decoderObj, romState, refState);
  auto problem = pressio::rom::lspg::createDefaultProblemSteady(
    appObj, decoderObj, romState, refState);

  const auto & currFomState = problem.currentFomState();

  // here, the fom state should be [8 8 8 ...]
  Eigen::VectorXd trueFomState(fomSize);
  trueFomState.setConstant(8.);
  if (!trueFomState.isApprox(currFomState)) checkStr = "FAILED";

  using solver_t = MyFakeSolver<rom_state_t,typename decoder_t::jacobian_type>;
  solver_t solver(fomSize, romSize, checkStr, currFomState);
  solver.solve(problem.getSystemRef(), romState);

  std::cout << checkStr <<  std::endl;
  return 0;
}
