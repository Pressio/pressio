
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_steady.hpp"

struct MyFom
{
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

  state_type createState() const{ return state_type(N_); }

  residual_type createResidual() const{ return residual_type(N_); }

  Eigen::MatrixXd createApplyJacobianResult(const Eigen::MatrixXd & B) const
  {
    Eigen::MatrixXd A(N_, B.cols());
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE(u.size()==r.size());
    EXPECT_TRUE(u.size()==N_);

    for (auto i=0; i<r.rows(); ++i){
     r(i) = u(i) + 1.;
    }
  }

  void applyJacobian(const state_type & state,
                     const Eigen::MatrixXd & B,
                     Eigen::MatrixXd & A) const
  {
    A = B;
    A.array() += 1.;
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

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, J);

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
          EXPECT_DOUBLE_EQ(J(i,j), start + (double)count++);
        }
      }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, J);

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
          EXPECT_DOUBLE_EQ(J(i,j), start + (double)count++);
        }
      }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

TEST(rom_lspg_steady, test1)
{

  /*
    - phi in R^{10,3}:
        phi[0,:]=0,1,2
        phi[1,:]=3,4,5
        phi[2,:]=6,7,8
        phi[3,:]=9,10,11
        phi[4,:]=12,13,14
        phi[5,:]=15,16,17
        phi[6,:]=18,19,20
        phi[7,:]=21,22,23

    - initial romState = [0,1,2]

    - fom residual R(y) always computes R[:] = y[:]+1
    - fom applyJac appJac(B) always returns B += 1
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  constexpr int N = 8;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
  int count = 0;
  for (int i=0; i<N; ++i){
    for (int j=0; j<3; ++j){
      phi(i,j) = (double) count++;
    }
  }
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = typename fom_t::state_type;
  auto space = pressio::rom::create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = pressio::rom::lspg::create_default_problem(space, fomSystem);

  FakeNonLinSolverSteady nonLinSolver(N);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
