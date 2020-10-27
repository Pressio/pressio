
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

  void residual(const state_type & state,
		residual_type & r) const
  {
    ++counter_;
    for (auto i =0; i<r.size(); ++i)
      r(i) = (double) i;
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     dense_matrix_type & A) const
  {
    ++counter1_;
    for (auto j=0; j<A.cols(); ++j)
      for (auto i=0; i<A.rows(); ++i)
	A(i,j) = (double) i;
  }
};

template<typename r_t, typename j_t>
struct MyFakeSolver
{
  int maskSize_ = {};
  int romSize_ = {};
  int callCounter_ = 0;
  r_t R_;
  j_t J_;
  std::string & checkString_;

  MyFakeSolver(int maskSize,
	       int romSize,
	       std::string & checkString)
    : maskSize_(maskSize),
      romSize_(romSize),
      R_(maskSize),
      J_(maskSize, romSize),
      checkString_(checkString)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++callCounter_;
    Eigen::VectorXd trueR(maskSize_);
    Eigen::MatrixXd trueJ(maskSize_, romSize_);

    sys.residual(state, R_);
    sys.jacobian(state, J_);
    std::cout << *R_.data() << std::endl;
    std::cout << *J_.data() << std::endl;

    if (R_.extent(0) != maskSize_) checkString_ = "FAILED";
    if (J_.extent(0) != maskSize_) checkString_ = "FAILED";
    if (J_.extent(1) != romSize_) checkString_ = "FAILED";

    if (callCounter_==1)
      {
	trueR << 0,1,3,4,6;
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

	trueJ << 0,0,0,1,1,1,3,3,3,4,4,4,6,6,6;
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }
  }
};

int main(int argc, char *argv[])
{
  /* verify correctness of sequence of calls to steady masked lspg */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;
  using masker_t	= Masker;
  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;
  constexpr int fomSize = 7;
  constexpr int romSize = 3;
  constexpr int maskSize= 5;

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

  // masker
  masker_t Masker({0,1,3,4,6});

  auto lspgProblem = pressio::rom::lspg::createMaskedProblemSteady
    (appObj, decoderObj, romState, refState, Masker);

  using solver_t = MyFakeSolver<rom_state_t,typename decoder_t::jacobian_type>;
  solver_t solver(maskSize, romSize, checkStr);
  solver.solve(lspgProblem.systemRef(), romState);

  std::cout << checkStr <<  std::endl;
  return 0;
}
