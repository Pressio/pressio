
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

constexpr int N = 7;

using FomStateType = Eigen::VectorXd;
using FomRhsType = Eigen::VectorXd;
using RomStateType = Eigen::VectorXd;
using RomRType = Eigen::VectorXd;
using RomJType = Eigen::MatrixXd;

struct MyFom
{
  using time_type = double;
  using state_type = FomStateType;
  using right_hand_side_type = FomRhsType;

  MyFom(){}


  right_hand_side_type createRightHandSide() const{
    right_hand_side_type r(N);
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

  Eigen::MatrixXd createResultOfJacobianActionOn(const Eigen::MatrixXd & A) const{
    return Eigen::MatrixXd(N, A.cols());
  }

  void applyJacobian(const state_type & s,
		     const Eigen::MatrixXd & A,
		     const time_type & evaltime,
		     Eigen::MatrixXd & result) const
  {
    result.col(0).setConstant(-2.);
    result.col(1).setConstant(-3.);
    result.col(2).setConstant(-4.);
  }

  Eigen::MatrixXd createResultOfMassMatrixActionOn(const Eigen::MatrixXd & operand) const{
    return Eigen::MatrixXd(N, operand.cols());
  }

  void applyMassMatrix(const Eigen::VectorXd & stateIn,
		       const Eigen::MatrixXd & operand,
		       double evalTime,
		       Eigen::MatrixXd & result) const
  {
    Eigen::MatrixXd M(N, N);
    for (std::size_t j=0; j<M.cols(); ++j){
      M.col(j) = stateIn;
      for (std::size_t i=0; i<M.rows(); ++i){
	M(i,j) += evalTime + (double) j;
      }
    }
    result = M * operand;
  }
};

struct NonLinSolver
{
  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    auto R = system.createResidual();
    auto J = system.createJacobian();

    // do solver iterator 1
    system.residualAndJacobian(state, R, J, true);
    // we fake an update to the state
    for (int i=0; i<state.size(); ++i){ state[i] += 1.; }
  }
};

TEST(rom_galerkin, test)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  // create fom
  using fom_t = MyFom;
  fom_t fomSystem;

  // create trial space
  using basis_t = Eigen::MatrixXd;
  basis_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);

  typename MyFom::state_type shift(N);
  auto space = pressio::rom::create_trial_column_subspace<RomStateType>(phi,
								 shift,
								 false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  const auto odeScheme = pressio::ode::StepScheme::BDF1;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_implicit_problem(odeScheme, space, fomSystem);

  using time_type = typename fom_t::time_type;
  const time_type dt = 2.;
  NonLinSolver solver;

  pressio::ode::advance_n_steps(problem, romState, time_type{0}, dt,
				::pressio::ode::StepCount(1), solver);
  std::cout << romState << std::endl;
  EXPECT_TRUE(false);

  pressio::log::finalize();
}
