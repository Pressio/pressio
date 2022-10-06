
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

struct Observer
{
  void operator()(pressio::ode::StepCount stepIn,
		  double /*time*/,
		  const Eigen::VectorXd & state) const
  {
    const auto step = stepIn.get();

    // EXPECT_TRUE(step<=2);

    // if (step==0){
    //   EXPECT_DOUBLE_EQ(state[0], 0.);
    //   EXPECT_DOUBLE_EQ(state[1], 1.);
    //   EXPECT_DOUBLE_EQ(state[2], 2.);
    // }
    // if (step==1){
    //   EXPECT_DOUBLE_EQ(state[0], 0.);
    //   EXPECT_DOUBLE_EQ(state[1], 51.);
    //   EXPECT_DOUBLE_EQ(state[2], 102.);
    // }
    // if (step==2){
    //   EXPECT_DOUBLE_EQ(state[0], 0.);
    //   EXPECT_DOUBLE_EQ(state[1], 2611.);
    //   EXPECT_DOUBLE_EQ(state[2], 5222.);
    // }
  }
};

struct MyFom
{
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

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
};

struct FomMassMatrixEvaluator
{
  int N_ = {};

  FomMassMatrixEvaluator(int N): N_(N){}

  Eigen::MatrixXd createApplyMassMatrixResult(const Eigen::MatrixXd & operand) const{
    return Eigen::MatrixXd(N_, operand.cols());
  }

  void operator()(const Eigen::VectorXd & stateIn,
		  const Eigen::MatrixXd & operand,
		  double evalTime,
		  Eigen::MatrixXd & result) const
  {
    Eigen::MatrixXd M(N_, N_);
    for (std::size_t j=0; j<M.cols(); ++j){
      M.col(j) = stateIn;
      for (std::size_t i=0; i<M.rows(); ++i){
	M(i,j) += evalTime + (double) j;
      }
    }
    result = M * operand;
  }
};

template<class MatrixType, class VecType>
struct FakeLinearSolver1
{
  void solve(const MatrixType & A, VecType & x, const VecType & b)
  {

    // we don't need to do anything meaningful here
    // we just need numbers. To make things simple
    // pretend that solving this system boils down
    // to doing x = Ab
    x = A*b;
    std::cout << x << std::endl;

    MatrixType goldA (3,3);
    goldA << 400., 800., 1200.,
      800., 1600., 2400., \
      1200., 2400., 3600;
    EXPECT_TRUE( goldA.isApprox(A) );

    VecType gold_b (3);
    gold_b << 70., 140., 210.;
    EXPECT_TRUE( gold_b.isApprox(b) );
  }
};

TEST(rom_galerkin_unsteady, test5)
{
  /* default galerkin explicit with euler forward

     M* delta = phi^T f

     where:
     - M is the "FOM" size = 5

     - phi[:,0]=1.; phi[:,1]=2.; phi[:,2]=3.;

     - M(y,time) such that M[:,0] = y+time; M[:,1] = y+time+1; ...

     M* =
	  [ 400  800 1200
	    800 1600 2400
	   1200 2400 3600]

     phi^T f = [70; 140; 210]
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  constexpr int N = 5;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  // create trial space
  using basis_t = Eigen::MatrixXd;
  basis_t phi(N, 3);
  phi.col(0).setConstant(1.);
  phi.col(1).setConstant(2.);
  phi.col(2).setConstant(3.);

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type shift(N);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=1.;
  romState[1]=2.;
  romState[2]=3.;

  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  FomMassMatrixEvaluator fomMMEval(N);
  auto problem = gal::create_unsteady_explicit_problem(odeScheme, space, fomSystem, fomMMEval);

  using v_t = typename decltype(problem)::state_type;
  using M_t = typename decltype(problem)::mass_matrix_type;
  FakeLinearSolver1<M_t, v_t> linSolver;

  using time_type = typename fom_t::time_type;
  const time_type dt = 2.;
  Observer obs;
  pressio::ode::advance_n_steps(problem, romState, time_type{0}, dt,
				::pressio::ode::StepCount(1), obs,
				linSolver);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 784001.);
  EXPECT_DOUBLE_EQ(romState[1], 1568002.);
  EXPECT_DOUBLE_EQ(romState[2], 2352003.);

  pressio::log::finalize();
}
