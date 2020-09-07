
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  mutable int stepsCount1_ = 0;
  mutable int stepsCount2_ = 0;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : romSize_{romSize}, jac_(fomSize, romSize)
  {
    //initialize the jacobian to be [1,2,3..]
    for (int i=0; i<fomSize; ++i){
      jac_(i,0) = 1; jac_(i,1) = 2; jac_(i,2) = 3; jac_(i,3) = 4;
    }
  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    ::pressio::containers::Vector<Eigen::VectorXd> & result) const
  {
    const auto & jacNativeObj = *jac_.data();
    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = jacNativeObj * romStateNativeObj;
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

  const jacobian_type & getReferenceToJacobian() const{
    return jac_;
  }
};

struct MyFakeApp
{
  const int N_ = {};
  std::string & checkStr_;
  mutable int solverCallCount_ = 0;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using discrete_time_residual_type = residual_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N, std::string & checkStr)
    : N_(N), checkStr_(checkStr){}

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
			    pressio::Norm normKind,
			    scalar_type & normR,
  			    Args && ... states) const
  {
    discreteTimeResidualImpl(step, time,dt,R,normKind,normR,
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
				pressio::Norm normKind,
				scalar_type & normR,
				const state_type & yn,
				const state_type & ynm1) const
  {
    ++solverCallCount_;

    std::cout << "step " << step << std::endl;
    std::cout << "yn " << std::endl;
    std::cout << yn << std::endl;
    std::cout << "ynm1 " << std::endl;
    std::cout << ynm1 << std::endl;

    if (solverCallCount_==1)
    {
      if (step!=1) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
	auto diff = std::abs(yn(i) - 0);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
	auto diff = std::abs(ynm1(i) - 0);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

    if (solverCallCount_==2)
    {
      if (step!=1) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
	auto diff = std::abs(yn(i) - 18.);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
	auto diff = std::abs(ynm1(i) - 0);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

    if (solverCallCount_==3)
    {
      if (step!=2) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
	auto diff = std::abs(yn(i) - 44.);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
	auto diff = std::abs(ynm1(i) - 44.);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

    if (solverCallCount_==4)
    {
      if (step!=2) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
	auto diff = std::abs(yn(i) - 78.);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
	auto diff = std::abs(ynm1(i) - 44.);
	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

  }
};

template<typename r_t, typename j_t>
struct MyFakeSolver
{
  int romSize_ = {};
  int callCounter_ = 0;
  r_t R_;
  j_t J_;
  Eigen::MatrixXd trueJ1;
  Eigen::MatrixXd trueJ2;
  Eigen::VectorXd trueR1;
  Eigen::VectorXd trueR2;
  std::string & checkString_;

  MyFakeSolver(int romSize, std::string & checkString)
    : romSize_(romSize),
      R_(romSize),
      J_(romSize, romSize),
      trueJ1(romSize_, romSize_),
      trueJ2(romSize_, romSize_),
      trueR1(romSize_),
      trueR2(romSize_),
      checkString_(checkString)
  {}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++callCounter_;

    // // here sys is a an "arbitrary" stepper object since for
    // // this test we use the discrete time API
    double norm{};
    for (auto k=0; k<2; ++k)
    {
      sys.residual(state, R_, ::pressio::Norm::L2, norm);
      sys.jacobian(state, J_);
    //   std::cout << *R_.data() << std::endl;

    //   // check that R and J are what they are supposed to be
    //   if (k==0 and callCounter_==1)
    //   {
    // 	for (int i=0; i<R_.extent(0); ++i){
    // 	  auto diff = std::abs(R_(i) - trueR1(i));
    // 	  //std::cout << diff << std::endl;
    // 	  if (diff > 1e-13 ) checkString_ = "FAILED";
    // 	}

    // 	for (int i=0; i<J_.extent(0); ++i){
    // 	  for (int j=0; j<J_.extent(1); ++j){
    // 	    auto diff = std::abs(J_(i,j) - trueJ1(i,j));
    // 	    if (diff > 1e-13 ) checkString_ = "FAILED";
    // 	  }
    // 	}
    //   }

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};

int main(int argc, char *argv[])
{
  /* In this test we solve:
       x_n+1 - n = dt * phi^T f( phi x_n+1)

     where x is vec of generalized coords.

     This test is meant to test Galerkin with the discrete time API
     works correctly when using a decoder where its Jacobian (J_d) changes.
     Obviously, we should use a solver to solve the system but we fake solver.
     To do this, we craft a test where:

     * dt = 0.1, we do steps so we go from t_0->t_1->t_2
     * the fake solver, which increments solution vector by 1.
     * R=[1 1 ... 1]^T always for simplicity

     --------
     STEP 1
     --------
     goes from t_0 to t_1:
     - start from x = [0 0 0 0]
       and J_d0 = [1 2 3 4;
		   1 2 3 4;
	           ...
		   1 2 3 4];

     When solver starts iter=1:
	J_d = J_d0 + 1
	x = [0 0 0 0], x_fom_n = [0 0... 0], and x_fom_nm1 = [0 0 ... 0]
     When solver starts iter=2:
	J_d = J_d0 + 2
	x = [1 1 1 1], x_fom_n = [18 18... 18], and x_fom_nm1 = [0 0 ... 0],

     at t_1 we should have: x = [2 2 2 ]

     --------
     STEP 2
     --------
     goes from t_1 to t_2:
     start from x = [2 2 2 2]
     and J_d0 = [3 4 5 6;
		 3 4 5 6;
	         ...
		 3 4 5 6];

     when solver starts iter=1:
	J_d = J_d0 + 1
	x = [2 2 2 2], x_fom_n = [44 44... 44], and x_fom_nm1 = [44 44 ... 44]
     when solver starts iter=2:
	J_d = J_d0 + 2
	x = [3 3 3 3], x_fom_n = [78 78... 78], and x_fom_nm1 = [44 44 ... 44],

   */


  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 10;
  constexpr int romSize = 4;

  // app object
  fom_t appObj(fomSize, checkStr);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  auto t0 = static_cast<scalar_t>(0);
  using ode_tag = pressio::ode::implicitmethods::Arbitrary;
  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using rom_jacobian_t = pressio::containers::Matrix<Eigen::MatrixXd>;

  using problem_t = pressio::rom::galerkin::composeDefaultProblem<
    ode_tag, fom_t, rom_state_t, rom_jacobian_t,
    decoder_t, stepper_order, stepper_n_states>::type;
  using stepper_t  = typename problem_t::stepper_t;
  problem_t Problem(appObj, refState, decoderObj, romState, t0);

  auto & stepperObj = Problem.getStepperRef();
  MyFakeSolver<rom_state_t, rom_jacobian_t> solver(romSize, checkStr);
  pressio::ode::advanceNSteps(stepperObj, romState, 0.0, 0.1, 2, solver);

  std::cout << *romState.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  return 0;
}
