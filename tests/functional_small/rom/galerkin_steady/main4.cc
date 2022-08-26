
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_steady.hpp"

struct MyFom
{
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  MyFom(int N,  std::vector<int> ind): N_(N), indices_to_corrupt_(ind){}

  state_type createState() const{ return state_type(N_); }
  residual_type createResidual() const{ return residual_type(N_); }

  Eigen::MatrixXd createApplyJacobianResult(const Eigen::MatrixXd & B) const{
    Eigen::MatrixXd A(N_, B.cols());
    return A;
  }

  void residual(const state_type & u, residual_type & r) const{
    EXPECT_TRUE(u.size()==r.size());
    EXPECT_TRUE(u.size()==N_);
    for (auto i=0; i<r.rows(); ++i){
     r(i) = u(i) + 1.;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     r(it) = -1114;
    }
  }

  void applyJacobian(const state_type & state,
                     const Eigen::MatrixXd & B,
                     Eigen::MatrixXd & A) const{
    A = B;
    for (int i=0; i<A.rows(); ++i){
      for (int j=0; j<A.cols(); ++j){
        A(i,j) += state(i);
      }
    }
    for (auto & it : indices_to_corrupt_){
      for (int j=0; j< A.cols(); ++j){
        A(it,j) = -4232;
      }
    }
  }
};

struct FakeNonLinSolverSteady
{
  int call_count_ = 0;
  int N_ = {};

  FakeNonLinSolverSteady(int N) : N_(N){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // mimic solver iterator 1
      system.residualAndJacobian(state, R, J, true);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;

      EXPECT_DOUBLE_EQ(R(0), 0.);
      EXPECT_DOUBLE_EQ(R(1), 48.);
      EXPECT_DOUBLE_EQ(R(2), 96.);

      Eigen::MatrixXd goldJ(3,3);
      goldJ.row(0) << 0., 0., 0.;
      goldJ.row(1) << 5.*8., 6.*8., 7*8.;
      goldJ.row(2) << 10.*8., 12.*8., 14*8.;
      EXPECT_TRUE(goldJ.isApprox(J));
    }

    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

    {
      // mimic solver iterator 2
      system.residualAndJacobian(state, R, J, true);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;

      EXPECT_DOUBLE_EQ(R(0), 0.);
      EXPECT_DOUBLE_EQ(R(1), 9.*8.);
      EXPECT_DOUBLE_EQ(R(2), 18.*8.);

      Eigen::MatrixXd goldJ(3,3);
      goldJ.row(0) << 0., 0., 0.;
      goldJ.row(1) << 8.*8, 9*8., 10*8.;
      goldJ.row(2) << 16.*8, 18*8., 20*8.;
      EXPECT_TRUE(goldJ.isApprox(J));

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

struct MyHypRedOperator
{
  using operator_type = Eigen::MatrixXd;
  operator_type matrix_;

public:
  using residual_operand_type = Eigen::VectorXd;
  using jacobian_action_operand_type = Eigen::MatrixXd;

  MyHypRedOperator(const operator_type & phiSampleMesh)
    : matrix_(phiSampleMesh){}

  template<class ResultType>
  void operator()(const residual_operand_type & operand, ResultType & result) const{
    result = matrix_.transpose() * operand;
  }

  template<class ResultType>
  void operator()(const jacobian_action_operand_type & operand, ResultType & result) const{
    result = matrix_.transpose() * operand;
  }
};

struct MyMaskerResidual
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = Eigen::VectorXd;
  using result_type = operand_type;

  MyMaskerResidual(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  result_type createApplyMaskResult(const operand_type & /*operand*/) const{
    return result_type(sample_indices_.size());
  }

  void operator()(const operand_type & operand, result_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

struct MyMaskerJacAction
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = Eigen::MatrixXd;
  using result_type = operand_type;

  MyMaskerJacAction(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

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

TEST(rom_galerkin_steady, test4)
{
  /* steady galerkin masked */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  constexpr int N = 13;
  /* corrupt indices are those that we mess up on purpose */\
  const std::vector<int> corrupt_indices = {1,3,4,5,9};\
  const std::vector<int> sample_indices = {0,2,6,7,8,10,11,12};\

  using fom_t = MyFom;
  fom_t fomSystem(N, corrupt_indices);

  /* phi must be valid to reconstruct fhe full fom */\
  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type shift(N);
  auto space = pressio::rom::create_trial_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  /* projector must be applicable to the *masked* operand*/\
  /* so we need to use only certain rows of phi*/\
  phi_t matForProj(N - corrupt_indices.size(), 3);\
  matForProj.col(0).setConstant(0.);
  matForProj.col(1).setConstant(1.);
  matForProj.col(2).setConstant(2.);
  MyHypRedOperator proj(matForProj);

  MyMaskerResidual  masker1(sample_indices);
  MyMaskerJacAction masker2(sample_indices);

  auto problem = pressio::rom::galerkin::create_steady_problem(space, fomSystem,
							       masker1, masker2,
							       proj);

  FakeNonLinSolverSteady nonLinSolver(N);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
