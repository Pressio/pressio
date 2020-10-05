
#include "pressio_rom.hpp"

constexpr double dt = 0.5;

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

  const jacobian_type & getReferenceToJacobian() const{
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
    ++callCounter_;

    // mimic the case where the velocity is computed
    // here the velocity, f, is always one
    Eigen::VectorXd f(N_);
    f.setConstant(1.);

    R = yn - ynm1 - dt*f;

    Eigen::VectorXd stateExpected(N_);
    if(callCounter_==1){
      stateExpected.setConstant(0.);
      if( ! stateExpected.isApprox(yn) ) checkStr_="FAILED";
    }
    if(callCounter_==2){
      stateExpected.setConstant(6.);
      if( ! stateExpected.isApprox(yn) ) checkStr_="FAILED";
    }
    if(callCounter_==3){
      stateExpected.setConstant(12.);
      if( ! stateExpected.isApprox(yn) ) checkStr_="FAILED";
    }
    if(callCounter_==4){
      stateExpected.setConstant(18.);
      if( ! stateExpected.isApprox(yn) ) checkStr_="FAILED";
    }
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
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(-dt);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 1. - dt*10;
      	  trueJ(i,1) = 2. - dt*20;
      	  trueJ(i,2) = 3. - dt*30;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==1 and k==1)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(6.-dt);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ.coeffRef(i,0) = 2. - dt*20.;
      	  trueJ.coeffRef(i,1) = 3. - dt*30.;
      	  trueJ.coeffRef(i,2) = 4. - dt*40.;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==2 and k==0)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(-dt);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 3. - dt*30;
      	  trueJ(i,1) = 4. - dt*40;
      	  trueJ(i,2) = 5. - dt*50;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }

      if (callCounter_==2 and k==1)
      {
      	Eigen::VectorXd trueR(fomSize_); trueR.setConstant(18.-12-dt);
      	if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      	for (auto i=0; i<trueJ.rows(); ++i){
      	  trueJ(i,0) = 4. - dt*40;
      	  trueJ(i,1) = 5. - dt*50;
      	  trueJ(i,2) = 6. - dt*60;
      	}
      	if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
      }


      for (auto i=0; i<state.extent(0); ++i) state(i) += 1.;

    }
  }
};


int main(int argc, char *argv[])
{
  /* Here we verify correctness of the varying decoder's jacobian
     for LSPG with discrete-time API.

     Let g(x) be the decoder, and let Jg be its jacobian.

     We solve: R() = dx/dt = Jg^T f( g(x) )
     using Euler forward where x is vec of generalized coords.

     We then have:
        R( ) = x_n+1 - x_n - dt * f( g(x_n+1) )
	J( ) = I - dt df/dx_fom Jg, where x_fom = g(x)

     To do this, we craft a test where:

     - dt = 0.5, we do 2 steps so from t_0->t_1->t_2

     - g(x) = [ 2 2 2 ]  x
                2 2 2
		...
	      [ 2 2 2 ]

     - f=[1 1 ... 1]^T always

     - df/dx_fom which is the app jacobian (Japp) is always
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]
        [ 1 2 1 2 1 2 1]

     --------
     STEP 1:
     --------
     * from t_0 to t_1
     * x_n = [0 0 0], x_nm1 = [0 0 0]

     When solver starts iter=1:
	Jg = [1 2 3]
	     |1 2 3|
	     |1 2 3|
	     ...
	     |1 2 3|

	x_n+1     = [0 0 0], and x_n = [0 0 0]
	x_fom_n+1 = [0 0... 0], and x_fom_n = [0 0 ... 0]

	R = -dt * f
	J = I*Jg - dt * Japp Jg

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

	x_n+1     = [1 1 1], and x_n = [0 0 0]
	x_fom_n+1 = [6 6... 6], and x_fom_n = [0 0 ... 0]

	R = y_fom_n+1 -dt * f
	J = I*Jg - dt * Japp Jg

	where Japp Jg =
        [ 1 2 ... ]   [2 3 4]
	| 1 2 ... |   |2 3 4|
	...         *
	| 1 2 ... |   |2 3 4|
	[ 1 2 ... ]   [2 3 4]


     --------
     STEP 2:
     --------
     * from t_1 to t_2

     When solver starts iter=1:
	Jg = [3 4 5]
	     |3 4 5|
	     |3 4 5|
	     ...
	     |3 4 5|

	x_n+1     = [2 2 2], and x_n = [2 2 2]
	x_fom_n+1 = [12 12... 12], and x_fom_n = [12 12 ... 12]

	R = -dt * f
	J = I*Jg - dt * Japp Jg

	where Japp Jg =
        [ 1 2 ... ]   [3 4 5]
	| 1 2 ... |   |3 4 5|
	...         *
	| 1 2 ... |   |3 4 5|
	[ 1 2 ... ]   [3 4 5]

     When solver starts iter=2:
	Jg = [4 5 6]
	     |4 5 6|
	     |4 5 6|
	     ...
	     |4 5 6|

	x_n+1     = [3 3 3], and x_n = [2 2 2]
	x_fom_n+1 = [18 18... 18], and x_fom_n = [12 12 ... 12]

	R = y_fom_n+1 -dt * f
	J = I*Jg - dt * Japp Jg

	where Japp Jg =
        [ 1 2 ... ]   [4 5 6]
	| 1 2 ... |   |4 5 6|
	...         *
	| 1 2 ... |   |4 5 6|
	[ 1 2 ... ]   [4 5 6]

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
  refState.setConstant(0.);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  using ode_tag = pressio::ode::implicitmethods::Arbitrary;
  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

  using problem_t  = pressio::rom::lspg::composeDefaultProblem<
    ode_tag, fom_t, rom_state_t, decoder_t, stepper_order, stepper_n_states>::type;
  problem_t prob(appObj, refState, decoderObj, romState);

  MyFakeSolver<rom_state_t, typename decoder_t::jacobian_type> solver(fomSize, romSize, checkStr);
  pressio::ode::advanceNSteps(prob.getStepperRef(), romState, 0.0, dt, 2, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
