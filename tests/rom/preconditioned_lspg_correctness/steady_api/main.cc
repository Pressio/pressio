
#include "pressio_rom.hpp"
#include "../helpers.hpp"

struct MyFakeApp
{
  int N_ = {};
  std::string & sentinel_;
  mutable int counter_ = {};
  mutable int counter1_ = {};

public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N, std::string & sentinel)
    : N_(N), sentinel_(sentinel){}

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
    ++counter1_;
    if (counter1_==1) A.setConstant(1);
    if (counter1_==2) A.setConstant(2);
  }

  void residual(const state_type & state,
		residual_type & r) const
  {
    ++counter_;
    if (counter_==1) r.setConstant(1);
    if (counter_==2) r.setConstant(2);
  }
};

template<typename r_t, typename j_t>
struct MyFakeSolver
{
  int fomSize_ = {};
  int romSize_ = {};
  int callCounter_ = 0;
  r_t R_;
  j_t J_;
  std::string & checkString_;

  MyFakeSolver(int fomSize, int romSize,
	       std::string & checkString)
    : fomSize_(fomSize),
      romSize_(romSize),
      R_(fomSize),
      J_(fomSize, romSize),
      checkString_(checkString)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++callCounter_;
    Eigen::MatrixXd trueJ(fomSize_, romSize_);

    for (auto k=0; k<2; ++k)
    {
      std::cout << "Solver call " << callCounter_ << " " << k << std::endl;

      sys.residual(state, R_);
      sys.jacobian(state, J_);

      std::cout << *state.data() << std::endl;
      std::cout << *R_.data() << std::endl;
      std::cout << *J_.data() << std::endl;

      if (callCounter_==1 and k==0)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(2.);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 1.+static_cast<double>(i);
	  trueJ(i,1) = 1.+static_cast<double>(i);
	  trueJ(i,2) = 1.+static_cast<double>(i);
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==1 and k==1)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(3.);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 2.+static_cast<double>(i);
	  trueJ(i,1) = 2.+static_cast<double>(i);
	  trueJ(i,2) = 2.+static_cast<double>(i);
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};

int main(int argc, char *argv[])
{
  /*  verify correctness of sequence of calls
      for steady preconditioned lspg.
      This test only checks that the preconditioner is called correctly.
   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;
  using precond_t	= Preconditioner;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 7;
  constexpr int romSize = 3;

  // app object
  fom_t appObj(fomSize, checkStr);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);
  refState.setConstant(0.0);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  // prec obj
  precond_t PrecObj;

  // using lspg_problem_type
  //   = typename pressio::rom::lspg::composePreconditionedProblem<
  //     fom_t, decoder_t, rom_state_t, precond_t>::type;
  // lspg_problem_type lspgProblem(appObj, decoderObj, romState, refState, PrecObj);
  auto lspgProblem = pressio::rom::lspg::createPreconditionedProblemSteady(
    appObj, decoderObj, romState, refState, PrecObj);
  using traits = ::pressio::rom::details::traits<decltype(lspgProblem)>;
  static_assert( std::is_same<typename traits::ud_ops_t, void>::value, "");

  using solver_t = MyFakeSolver<rom_state_t,typename decoder_t::jacobian_type>;
  solver_t solver(fomSize, romSize, checkStr);
  solver.solve(lspgProblem.systemRef(), romState);

  std::cout << checkStr <<  std::endl;
  return 0;
}
