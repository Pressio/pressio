
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_steady.hpp"

struct MyFom
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  int nSample_  = {};
  int nStencil_ = {};
  const std::vector<int> indices_ = {0,2,4,6,8,10,12,14};

  MyFom(int nSample, int nStencil)
    : nSample_(nSample),
      nStencil_(nStencil){}

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

    for (std::size_t i=0; i<indices_.size(); ++i){
      r(i) = u(indices_[i]) + 1.;
    }

    if (Ain.value()){
      auto & A = *Ain.value();
      for (std::size_t i=0; i<indices_.size(); ++i){
	for (int j=0; j< A.cols(); ++j){
	  A(i,j) = B(indices_[i], j);
	  A(i,j) += 1.;
	}
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
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);
    using Jo_t = std::optional<decltype(J) *>;

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, Jo_t(&J));

      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R(0), 0*0.+1*1.+2*2.+1.);
      EXPECT_DOUBLE_EQ(R(1), 3*0.+4*1.+5*2.+1.);
      EXPECT_DOUBLE_EQ(R(2), 6*0.+7*1.+8*2.+1.);
      EXPECT_DOUBLE_EQ(R(3), 9*0.+10*1.+11*2.+1.);
      EXPECT_DOUBLE_EQ(R(4), 12*0.+13*1.+14*2.+1.);
      EXPECT_DOUBLE_EQ(R(5), 15*0.+16*1.+17*2.+1.);
      EXPECT_DOUBLE_EQ(R(6), 18*0.+19*1.+20*2.+1.);
      EXPECT_DOUBLE_EQ(R(7), 21*0.+22*1.+23*2.+1.);

      double start = 1;
      int count = 0;
      for (int i=0; i<N_; ++i){
        for (int j=0; j<3; ++j){
          EXPECT_NEAR(J(i,j), start + (double)count++, 1e-15);
        }
      }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, Jo_t(&J));

      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R(0), 0*1.+1*2.+2*3.+1.);
      EXPECT_DOUBLE_EQ(R(1), 3*1.+4*2.+5*3.+1.);
      EXPECT_DOUBLE_EQ(R(2), 6*1.+7*2.+8*3.+1.);
      EXPECT_DOUBLE_EQ(R(3), 9*1.+10*2.+11*3.+1.);
      EXPECT_DOUBLE_EQ(R(4), 12*1.+13*2.+14*3.+1.);
      EXPECT_DOUBLE_EQ(R(5), 15*1.+16*2.+17*3.+1.);
      EXPECT_DOUBLE_EQ(R(6), 18*1.+19*2.+20*3.+1.);
      EXPECT_DOUBLE_EQ(R(7), 21*1.+22*2.+23*3.+1.);

      start = 1;
      count = 0;
      for (int i=0; i<N_; ++i){
        for (int j=0; j<3; ++j){
          EXPECT_NEAR(J(i,j), start + (double)count++, 1e-15);
        }
      }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

    }// end if call_count == 1
  }// end solve
};

TEST(rom_lspg_steady, test3)
{
  /*
    steady hyper-reduced lspg
   */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  const int nStencil = 15;
  const int nSample  = 8;

  using fom_t = MyFom;
  fom_t fomSystem(nSample, nStencil);

  Eigen::MatrixXd phi(nStencil, 3);
  int count = 0;
  for (int i=0; i<nStencil; ++i){
    for (int j=0; j<3; ++j){
      if (i % 2 == 0){
        phi(i,j) = (double) count++;
      }
      else{
        phi(i,j) = (double) -1;
      }
    }
  }
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type dummyFomState(nStencil);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, dummyFomState, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = pressio::rom::lspg::create_steady_problem(space, fomSystem);

  FakeNonLinSolverSteady nonLinSolver(nSample);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
