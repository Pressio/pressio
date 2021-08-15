
#include "pressio_rom_lspg.hpp"
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
    // fake modification
    if(callCounter1_==1) A.setConstant(1.);
    if(callCounter1_==2) A.setConstant(2.);
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
    // fake modification
    if(callCounter_==1) R.setConstant(1.);
    if(callCounter_==2) R.setConstant(2.);
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

    std::cout << "Solver call " << callCounter_ << " " << std::endl;

    sys.residual(state, R_);
    sys.jacobian(state, J_);

    std::cout << *state.data() << std::endl;
    std::cout << *R_.data() << std::endl;
    std::cout << *J_.data() << std::endl;

    if (callCounter_==1)
    {
      Eigen::VectorXd trueR(fomSize_); trueR.setConstant(2.);
      if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      for (auto i=0; i<trueJ.rows(); ++i){
      	trueJ(i,0) = static_cast<double>(i) + 1.;
      	trueJ(i,1) = static_cast<double>(i) + 1.;
      	trueJ(i,2) = static_cast<double>(i) + 1.;
      }
      if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
    }

    if (callCounter_==2)
    {
      Eigen::VectorXd trueR(fomSize_); trueR.setConstant(3.);
      if (!trueR.isApprox(*R_.data())) checkString_ = "FAILED";

      for (auto i=0; i<trueJ.rows(); ++i){
      	trueJ(i,0) = static_cast<double>(i) + 2.;
      	trueJ(i,1) = static_cast<double>(i) + 2.;
      	trueJ(i,2) = static_cast<double>(i) + 2.;
      }
      if (!trueJ.isApprox(*J_.data())) checkString_ = "FAILED";
    }
  }
};


int main(int argc, char *argv[])
{
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

  // using odetag = pressio::ode::implicitmethods::Arbitrary;
  // using stepper_order    = ::pressio::ode::StepperOrder<1>;
  // using stepper_n_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
  // using problem_t =
  //   pressio::rom::lspg::composePreconditionedDefaultProblem<
  //     ode_tag, fom_t, decoder_t, rom_state_t, precond_t,
  //   stepper_order, stepper_n_states>::type;
  // problem_t problem(appObj, decoderObj, romState, refState, PrecObj);
  auto problem = pressio::rom::lspg::createPreconditionedDefaultProblemUnsteady<1,2>(
    appObj, decoderObj, romState, refState, PrecObj);

  MyFakeSolver<rom_state_t, typename decoder_t::jacobian_type> solver(fomSize, romSize, checkStr);
  pressio::rom::lspg::solveNSequentialMinimizations(problem, romState, 0.0, dt, 2, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
