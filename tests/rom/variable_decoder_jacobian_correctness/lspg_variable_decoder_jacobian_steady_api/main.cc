
#include "pressio_rom_lspg.hpp"

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  mutable int applyMappingCount_ = 0;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : romSize_{romSize}, jac_(fomSize, romSize)
  {
    for (int i=0; i<fomSize; ++i){
      jac_(i,0) = 0; jac_(i,1) = 1; jac_(i,2) = 2;
    }
  }

  const jacobian_type & jacobianCRef() const{
    return jac_;
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) const
  {
    std::cout << "UPDATE\n";
    for (int i=0; i<jac_.extent(0); ++i){
      for (int j=0; j<jac_.extent(1); ++j){
	jac_(i,j) += 1;
      }
    }
  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    ::pressio::containers::Vector<Eigen::VectorXd> & result) const
  {
    ++applyMappingCount_;

    Eigen::MatrixXd A(result.extent(0), romSize_);
    A.setConstant(2);

    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = A * romStateNativeObj;
  }
};

struct MyFakeApp
{
  int N_ = {};
  std::string & sentinel_;
  mutable int counter_ = {};

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
    Eigen::MatrixXd J(N_, N_);
    for (auto i=0; i<N_; ++i){
      for (auto j=0; j<N_; j+=2) J(i,j) = 1.;
      for (auto j=1; j<N_; j+=2) J(i,j) = 2.;
    }
    A = J*B;
  }

  void residual(const state_type & state,
		residual_type & r) const
  {
    ++counter_;
    if(counter_==1) r.setConstant(1.);
    if(counter_==2) r.setConstant(2.);
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
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(1.);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 10;
      	  trueJ(i,1) = 20;
      	  trueJ(i,2) = 30;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==1 and k==1)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(2.);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 20;
      	  trueJ(i,1) = 30;
      	  trueJ(i,2) = 40;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};

int main(int argc, char *argv[])
{
  /* verify correctness of the varying decoder's jacobian for steady LSPG

     Let g(x) be the decoder, and let Jg be its jacobian.

     We solve: f( g(x) ) = 0

     We then have:
        R( ) = f( g(x_n+1) )
	J( ) = df/dx_fom * Jg, where x_fom = g(x)

     To do this, we craft a test where:

     - g(x) = [ 2 2 2 ]  x
                2 2 2
		...
	      [ 2 2 2 ]

     - f=[1 1 ... 1]^T at first call, then it is [2 2 ... 2]^T

     - df/dx_fom which is the app jacobian (Japp) is always
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]

     * x = [0 0 0]

     When solver starts iter=1:
	Jg = [1 2 3]
	     |1 2 3|
	     |1 2 3|
	     ...
	     |1 2 3|

	x = [0 0 0]
	x_fom = [0 0... 0]

	R = [1 1 ... 1]
	J = Japp Jg

	where Japp Jg =
        [ 1 2 ... ]   [1 2 3]
	| 1 2 ... |   |1 2 3|
	...         *
	| 1 2 ... |   |1 2 3|
	[ 1 2 ... ]   [1 2 3]

     When solver starts iter=2:
	Jg = [2 3 4]
	     |2 3 4|
	     ...
	     |2 3 4]

	x = [1 1 1]
	x_fom = [6 6... 6]

	R = [2 2 ... 2]
	J = Japp Jg

	where Japp Jg =
        [ 1 2 ... ]   [2 3 4]
	| 1 2 ... |   |2 3 4|
	...         *
	| 1 2 ... |   |2 3 4|
	[ 1 2 ... ]   [2 3 4]

   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;

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

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  // using lspg_problem_type = typename pressio::rom::lspg::composeDefaultProblem<
  //     fom_t, decoder_t, rom_state_t>::type;
  // lspg_problem_type lspgProblem(appObj, decoderObj, romState, refState);
  auto lspgProblem = pressio::rom::lspg::createDefaultProblemSteady(
    appObj, decoderObj, romState, refState);

  using solver_t = MyFakeSolver<rom_state_t,typename decoder_t::jacobian_type>;
  solver_t solver(fomSize, romSize, checkStr);
  solver.solve(lspgProblem.systemRef(), romState);

  std::cout << checkStr <<  std::endl;
  return 0;
}
