
#include "pressio_rom.hpp"
#include "../helpers.hpp"

constexpr double dt = 0.5;

struct MyFakeApp
{
  const int N_ = {};
  std::string & checkStr_;
  mutable int callCounter_ = 0;
  mutable int callCounter1_ = 0;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using discrete_time_residual_type = residual_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N, std::string & checkStr)
    : N_(N), checkStr_(checkStr){}

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
    ++callCounter1_;

    dense_matrix_type JB(N_,B.cols());
    for (auto j=0; j<JB.cols(); ++j){
      for (auto i=0; i<JB.rows(); ++i)
      {
	JB(i,j) = (double) i;
      }
    }
    if (callCounter1_==2) JB.array() += 1.;

    A = B - dt*JB;
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

private:
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    Eigen::VectorXd f(R.size());
    for (auto i =0; i<f.size(); ++i)
      f(i) = (double) i;

    for (auto i =0; i<R.size(); ++i)
      R(i) = yn(i) - ynm1(i) - dt*f(i);
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

    for (auto k=0; k<2; ++k)
    {
      std::cout << "state" << std::endl;
      std::cout << *state.data() << std::endl;
      sys.residual(state, R_);
      sys.jacobian(state, J_);
      std::cout << "residual" << std::endl;
      std::cout << *R_.data() << std::endl;
      std::cout << "jacobian" << std::endl;
      std::cout << *J_.data() << std::endl;

      if (R_.extent(0) != maskSize_) checkString_ = "FAILED";
      if (J_.extent(0) != maskSize_) checkString_ = "FAILED";
      if (J_.extent(1) != romSize_) checkString_ = "FAILED";

      if (callCounter_==1 and k==0)
      {
	trueR << 0, -0.5, -1.5, -2, -3;
	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

	trueJ <<
	  1-dt*0, 1-dt*0, 1-dt*0,
	  1-dt*1, 1-dt*1, 1-dt*1,
	  1-dt*3, 1-dt*3, 1-dt*3,
	  1-dt*4, 1-dt*4, 1-dt*4,
	  1-dt*6, 1-dt*6, 1-dt*6;
	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==2 and k==1)
      {
	trueR << 3, 3-0.5, 3-1.5, 3-2, 0.;
	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

	trueJ <<
	  1-dt*1, 1-dt*1, 1-dt*1,
	  1-dt*2, 1-dt*2, 1-dt*2,
	  1-dt*4, 1-dt*4, 1-dt*4,
	  1-dt*5, 1-dt*5, 1-dt*5,
	  1-dt*7, 1-dt*7, 1-dt*7;
	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      // fake solution
      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};


int main(int argc, char *argv[])
{
  /*
    ----------------------
    first call to solver:
    ----------------------
    R = y_n - y_n-1 -dt*f
      = -dt*f
      = [0, -0.5, -1.5, -2, -3.]

    J*phi = I*phi - dt*df/dy*phi
      where phi = [1 1 1;
		   ...
			]
     df/dy*phi = [ 0 0 0
                   1 1 1
		   ....
     only picked at sample mesh points

    ----------------------
    second call to solver:
    ----------------------
    R = y_n - y_n-1 -dt*f
      = [3 ... 3] -dt*f
      = [3-dt*0, 3-dt*1, 3-dt*3, 3-dt*4, 3-dt*6]

    J = I*phi - dt*df/dy*phi
      where phi = [1 1 1;
		   ...
			]
     df/dy*phi = [ 1 1 1
		   2 2 2
		   ....
     only picked at sample mesh points
   */

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

  auto problem = pressio::rom::lspg::createMaskedProblemUnsteady<1,2>
    (appObj, decoderObj, romState, refState, Masker);
  MyFakeSolver<rom_state_t, typename decoder_t::jacobian_type> solver
    (maskSize, romSize, checkStr);

  pressio::ode::advanceNSteps(problem.stepperRef(), romState, 0.0, dt, 1, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
