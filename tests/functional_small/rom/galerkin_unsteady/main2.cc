
#include <gtest/gtest.h>
#include "pressio/rom_tmp.hpp"

struct Observer
{
  void operator()(pressio::ode::StepCount stepIn,
		  double time,
		  const Eigen::VectorXd & state) const
  {
    const auto step = stepIn.get();

    EXPECT_TRUE(step<=2);

    if (step==0){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 1.);
      EXPECT_DOUBLE_EQ(state[2], 2.);
    }
    if (step==1){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 51.);
      EXPECT_DOUBLE_EQ(state[2], 102.);
    }
    if (step==2){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 2611.);
      EXPECT_DOUBLE_EQ(state[2], 5222.);
    }
  }
};

struct MyFom
{
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

  state_type createState() const{
    state_type s(N_);
    s.setConstant(0);
    return s;
  }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type r(N_);
    r.setConstant(0);
    return r;
  }

  void rightHandSide(const state_type & u,
		     const time_type evalTime,
		     right_hand_side_type & f) const
  {
    for (decltype(f.rows()) i=0; i<f.rows(); ++i){
      f(i) = u(i) + evalTime;
    }
  }

  Eigen::MatrixXd createApplyMassMatrixResult(const Eigen::MatrixXd & A) const{
    return Eigen::MatrixXd(N_, A.cols());
  }

  void applyMassMatrix(const state_type & s,
		       const Eigen::MatrixXd & A,
		       const time_type & evaltime,
		       Eigen::MatrixXd & result) const
  {
    std::cout << "gigi\n";
  }
};

struct LinearSolver1{
  void solve(const Eigen::MatrixXd & A,
       Eigen::VectorXd & x,
       const Eigen::VectorXd & b)
  {}
};

TEST(rom_galerkin, test2)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  // create fom
  constexpr int N = 10;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  // create trial space
  using basis_t = Eigen::MatrixXd;
  basis_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);

  using latent_state_type = Eigen::VectorXd;
  using full_state_type = typename fom_t::state_type;
  auto trialSpace = pressio::rom::create_linear_subspace<latent_state_type, full_state_type>(phi);

  auto romState = trialSpace.createLatentState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_default_explicit_problem(odeScheme, trialSpace, fomSystem);

  using time_type = typename fom_t::time_type;
  const time_type dt = 1.;
  LinearSolver1 ls;
  pressio::ode::advance_n_steps(problem, romState, time_type{0}, dt,
        ::pressio::ode::StepCount(2), ls);
  std::cout << romState << std::endl;

  // fix all checks
  EXPECT_TRUE(false);

  pressio::log::finalize();
}