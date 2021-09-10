
#include "pressio_rom_galerkin.hpp"

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;
  using fom_state_type = pressio::containers::Vector<Eigen::VectorXd>;

private:
  const int romSize_ = {};
  jacobian_type jac_;

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
  void applyMapping(const rom_state_type & romState, fom_state_type & result) const
  {
    Eigen::MatrixXd A(result.extent(0), romSize_);
    A.setConstant(2);

    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = A * romStateNativeObj;
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &)
  {
    std::cout << "UPDATE\n";
    for (int i=0; i<jac_.extent(0); ++i){
      for (int j=0; j<jac_.extent(1); ++j){
	jac_(i,j) += 1;
      }
    }
  }

  const jacobian_type & jacobianCRef() const{
    return jac_;
  }
};

struct MyFakeApp
{
  const int N_ = {};
  std::string & checkStr_;
  mutable int callCounter_ = 0;

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
    ++callCounter_;

    R.setConstant(1);

    std::cout << "step " << step << std::endl;
    std::cout << "yn " << std::endl;
    std::cout << yn << std::endl;
    std::cout << "ynm1 " << std::endl;
    std::cout << ynm1 << std::endl;

    if (callCounter_==1)
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

    if (callCounter_==2)
    {
      if (step!=1) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(yn(i) - 8.);
    	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(ynm1(i) - 0);
    	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

    if (callCounter_==3)
    {
      if (step!=2) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(yn(i) - 16.);
    	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(ynm1(i) - 16.);
    	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
    }

    if (callCounter_==4)
    {
      if (step!=2) checkStr_ = "FAILED";

      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(yn(i) - 24.);
    	if (diff > 1e-13 ) checkStr_ = "FAILED";
      }
      for (int i=0; i<N_; ++i){
    	auto diff = std::abs(ynm1(i) - 16.);
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
    for (auto k=0; k<2; ++k)
    {
      sys.residual(state, R_);
      sys.jacobian(state, J_);

      if (callCounter_==1 and k==0){
	Eigen::VectorXd trueR(romSize_);
	trueR << 20., 30., 40., 50.;
	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";
      }
      if (callCounter_==1 and k==1){
	Eigen::VectorXd trueR(romSize_);
	trueR << 30., 40., 50., 60.;
	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==2 and k==0){
	Eigen::VectorXd trueR(romSize_);
	trueR << 40., 50., 60, 70;
	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";
      }
      if (callCounter_==2 and k==1){
      	Eigen::VectorXd trueR(romSize_);
      	trueR << 50., 60., 70., 80.;
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";
      }

      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;
    }
  }
};

int main(int argc, char *argv[])
{
  /* Here we verify correctness of the varying decoder's jacobian
     for galerkin with discrete-time API.

     Let g(x) be the decoder, and let Jg be its jacobian.

     This test is meant to test Galerkin with the discrete time API
     works correctly when using a decoder where its Jacobian (Jg) changes.
     Obviously, we should use a solver to solve the system but we fake solver.
     To do this, we craft a test where:

     * dt = 0.1, we do steps so we go from t_0->t_1->t_2
     * the fake solver, which increments solution vector by 1.
     * R=[1 1 ... 1]^T always for simplicity

     - g(x) = [ 2 2 2 2 ]  x
                2 2 2 2
		...
	      [ 2 2 2 2 ]

     --------
     STEP 1
     --------
     - goes from t_0 to t_1:
     - start from x = [0 0 0 0]
       and Jg0 = [1 2 3 4;
		  1 2 3 4;
	          ...
		  1 2 3 4];

     When solver starts iter=1:
	Jg = Jg0 + 1
	x = [0 0 0 0], x_fom_n = [0 0... 0], and x_fom_nm1 = [0 0 ... 0]
	R[:] = 1
	Jg^T R = [ 2 2 ... 2]  R  = [ 20, 30, 40, 50]^T
		 [ 3 3 ... 3]
		 [ 4 4 ... 4]
		 [ 5 5 ... 5]

     When solver starts iter=2:
	Jg = Jg0 + 2
	x = [1 1 1 1], x_fom_n = [8 8... 8], and x_fom_nm1 = [0 0 ... 0],
	R[:] = 1
	Jg^T R = [ 3 3 ... 3]  R  = [ 30, 40, 50, 60]^T
		 [ 4 4 ... 4]
		 [ 5 5 ... 5]
		 [ 6 6 ... 6]

     at t_1 we should have: x = [2 2 2 ]

     --------
     STEP 2
     --------
     goes from t_1 to t_2:
     start from x = [2 2 2 2]
     and Jg0 = [3 4 5 6;
		3 4 5 6;
	        ...
		3 4 5 6];

     when solver starts iter=1:
	Jg = Jg0 + 1
	x = [2 2 2 2], x_fom_n = [16 16... 16], and x_fom_nm1 = [16 16 ... 16]
	Jg^T R = [ 4 4 ... 4]  R  = [ 40, 50, 60, 70]^T
		 [ 5 5 ... 5]
		 [ 6 6 ... 6]
		 [ 7 7 ... 7]

     when solver starts iter=2:
	Jg = Jg0 + 2
	x = [3 3 3 3], x_fom_n = [24 24... 24], and x_fom_nm1 = [16 16 ... 16],
	Jg^T R = [ 5 5 ... 5]  R  = [ 50, 60, 70, 80]^T
		 [ 6 6 ... 6]
		 [ 7 6 ... 7]
		 [ 8 7 ... 8]

   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
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
  refState.setConstant(0.0);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  using rom_jacobian_t = pressio::containers::DenseMatrix<Eigen::MatrixXd>;
  // using ode_tag = pressio::ode::ImplicitArbitrary;
  // using stepper_order    = ::pressio::ode::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
  // using problem_t = pressio::rom::galerkin::composeDefaultProblem<
  //   ode_tag, fom_t, decoder_t, rom_state_t, rom_jacobian_t,
  //   stepper_order, stepper_n_states>::type;
  // problem_t Problem(appObj, decoderObj, romState, refState);
  auto Problem =
    pressio::rom::galerkin::createDefaultProblem<rom_jacobian_t,1,2>(appObj, decoderObj, romState, refState);

  MyFakeSolver<rom_state_t, rom_jacobian_t> solver(romSize, checkStr);
  pressio::rom::galerkin::solveNSteps(Problem, romState, 0.0, 0.1, 2, solver);

  //std::cout << *romState.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  return 0;
}
