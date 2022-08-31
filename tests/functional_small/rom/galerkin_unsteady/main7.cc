
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

struct MyFom
{
  using time_type       = double;
  using state_type        = Eigen::VectorXd;
  using right_hand_side_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  MyFom(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    return OperandType(B.rows(), B.cols());
  }

  // computes: A = Jac B
  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     const time_type & time,
                     OperandType & A) const
  {
    A = B;
    A.array() += time;
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
    }
  }

  void rightHandSide(const state_type & u,
		     const time_type timeIn,
		     right_hand_side_type & f) const
  {
    for (auto i=0; i<f.rows(); ++i){
     f(i) = u(i) + timeIn;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
    }
  }
};

struct MyRhsMasker
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = typename MyFom::right_hand_side_type;
  using result_type = Eigen::VectorXd;

  MyRhsMasker(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  result_type createApplyMaskResult(const operand_type & operand) const{
    return result_type(sample_indices_.size());
  }

  void operator()(const operand_type & operand, result_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

struct MyJaMasker
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = Eigen::MatrixXd;
  using result_type = operand_type;

  MyJaMasker(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  result_type createApplyMaskResult(const operand_type & operand) const{
    return result_type(sample_indices_.size(), operand.cols());
  }

  void operator()(const operand_type & operand, result_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.cols(); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

class HypRedOperator
{
  using matrix_type = Eigen::MatrixXd;
  matrix_type matrix_;

public:
  using time_type = double;
  using right_hand_side_operand_type = typename MyRhsMasker::result_type;
  using jacobian_action_operand_type = typename MyJaMasker::result_type;

  HypRedOperator(const matrix_type & phi) : matrix_(phi){}

  void operator()(const right_hand_side_operand_type & operand,
		  time_type timein,
		  Eigen::VectorXd & result) const
  {
    result = matrix_.transpose() * operand;
  }

  void operator()(const jacobian_action_operand_type & operand,
		  time_type timein,
		  Eigen::MatrixXd & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

struct NonLinSolver
{
  int call_count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, J, true);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 0.);
      EXPECT_DOUBLE_EQ(R[1], -140.);
      EXPECT_DOUBLE_EQ(R[2], -280.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, J, true);
      EXPECT_DOUBLE_EQ(R[0], 1.);
      EXPECT_DOUBLE_EQ(R[1], -199.);
      EXPECT_DOUBLE_EQ(R[2], -399.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }

    //
    // call_count == 2
    //
    if(call_count_==2)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, J, true);
      EXPECT_DOUBLE_EQ(R[0],    0.);
      EXPECT_DOUBLE_EQ(R[1], -300.);
      EXPECT_DOUBLE_EQ(R[2], -600.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, J, true);
      EXPECT_DOUBLE_EQ(R[0],    1.);
      EXPECT_DOUBLE_EQ(R[1], -359.);
      EXPECT_DOUBLE_EQ(R[2], -719.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

TEST(rom_galerkin, test7)
{
  // test for masked implicit galerkin using BDF1
  // all numbers have been computed manually

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  constexpr int nFull = 20;
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14,16,18};
  const int nMasked = sample_indices.size();
  const std::vector<int> corrupt_indices = {1,7,13,19};

  MyFom fomSystem(nFull, corrupt_indices);

  using phi_t = Eigen::MatrixXd;
  phi_t phiFull(nFull, 3);
  phiFull.col(0).setConstant(0.);
  phiFull.col(1).setConstant(1.);
  phiFull.col(2).setConstant(2.);
  /* corrupting the entries in phiFull is the way we use to
     ensure the masking works properly */
  phiFull.row(1).setConstant(-111.);
  phiFull.row(7).setConstant(111.);
  phiFull.row(11).setConstant(423.);
  phiFull.row(17).setConstant(-21.);

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift(nFull);
  auto space = pressio::rom::create_trial_subspace<reduced_state_type>(phiFull,
								       shift,
								       false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  /* projector must be applicable to the *masked* operand*/
  /* so we need to use only certain rows of phi*/
  phi_t matForProj(nMasked, 3);
  for (int i = 0; i < nMasked; ++i){
    matForProj(i, 0) = phiFull(sample_indices[i],0);
    matForProj(i, 1) = phiFull(sample_indices[i],1);
    matForProj(i, 2) = phiFull(sample_indices[i],2);
  }\
  HypRedOperator hrOp(matForProj);

  MyRhsMasker rhsM(sample_indices);
  MyJaMasker jaM(sample_indices);

  const auto odeScheme = pressio::ode::StepScheme::BDF1;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_implicit_problem(odeScheme, space, fomSystem,
						       rhsM, jaM, hrOp);

  const double dt = 2.;
  NonLinSolver solver;
  pressio::ode::advance_n_steps(problem, romState, 0., dt,
				::pressio::ode::StepCount(2), solver);
  std::cout << romState << std::endl;
  EXPECT_TRUE(romState[0] == 4.);
  EXPECT_TRUE(romState[1] == 5.);
  EXPECT_TRUE(romState[2] == 6.);

  pressio::log::finalize();
}
