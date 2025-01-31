
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_steady.hpp"

namespace{

struct MyFom
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  int nSample_  = {};
  int nStencil_ = {};
  const std::vector<int> validStateIndices_{};

  MyFom(int nSample, int nStencil, const std::vector<int> & rowInd)
    : nSample_(nSample),
      nStencil_(nStencil),
      validStateIndices_(rowInd){}

  residual_type createResidual() const{ return residual_type(nSample_); }

  Eigen::MatrixXd createResultOfJacobianActionOn(const Eigen::MatrixXd & B) const{
    Eigen::MatrixXd A(nSample_, B.cols());
    return A;
  }

  void residualAndJacobianAction(const state_type & u,
				 residual_type & r,
				 const Eigen::MatrixXd & B,
				 std::optional<Eigen::MatrixXd *> Ain) const
  {
    EXPECT_TRUE(u.size()!=r.size());
    EXPECT_TRUE(u.size()==nStencil_);
    EXPECT_TRUE(r.size()==nSample_);

    for (std::size_t i=0; i<validStateIndices_.size(); ++i){
      r(i) = u(validStateIndices_[i]) + 1.;
    }

    if (Ain){
      auto & A = *Ain.value();
      for (std::size_t i=0; i<validStateIndices_.size(); ++i){
        for (int j=0; j< A.cols(); ++j){
          A(i,j) = B(validStateIndices_[i], j);
          A(i,j) += u(validStateIndices_[i]);
	      }
      }
    }
  }
};

struct FakeNonLinSolverSteady
{
  int call_count_ = 0;

  FakeNonLinSolverSteady(){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    if(call_count_==1)
    {
      // mimic solver iterator 1
      system.residualAndJacobian(state, R, &J);
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
      system.residualAndJacobian(state, R, &J);
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

class MyHypRedOperator
{
  using operator_type = Eigen::MatrixXd;
  operator_type matrix_;

public:
  MyHypRedOperator(const operator_type & phiSampleMesh)
    : matrix_(phiSampleMesh){}

  template<class T1, class T2>
  void operator()(const Eigen::MatrixBase<T1> & operand,
		  Eigen::MatrixBase<T2> & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

}

TEST(rom_galerkin_steady, hyperreduced)
{
  /*
    this test is numerically equivalent to main1
    except that we pretend here to do a hyper-reduced problem
   */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug, pressiolog::LogTo::console);

  const int nStencil = 13;
  const std::vector<int> validStateIndices = {2,3,4,5,6,10,11,12};
  const int nSample  = validStateIndices.size();

  using fom_t = MyFom;
  fom_t fomSystem(nSample, nStencil, validStateIndices);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(nStencil, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);
  phi.row(0).setConstant(-111.);
  phi.row(1).setConstant(-111.);
  phi.row(7).setConstant(111.);
  phi.row(8).setConstant(423.);
  phi.row(9).setConstant(-21.);
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type shift(nStencil);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  phi_t matForProj(nSample, 3);\
  matForProj.col(0).setConstant(0.);
  matForProj.col(1).setConstant(1.);
  matForProj.col(2).setConstant(2.);
  MyHypRedOperator hypRedOp(matForProj);
  auto problem = pressio::rom::galerkin::create_steady_problem(space, fomSystem, hypRedOp);

  FakeNonLinSolverSteady nonLinSolver;
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  PRESSIOLOG_FINALIZE();
}
