
#include <gtest/gtest.h>

#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"

struct JacobianAction{
  Eigen::VectorXd u_;
  const Eigen::MatrixXd * operand_;
  JacobianAction(int N, const Eigen::MatrixXd & B) :
    u_(N), operand_(&B)
  {
    u_.setConstant(0);
  }

  void store(const Eigen::VectorXd & u,
	     const Eigen::MatrixXd & newB)
  {
    u_ = u;
    operand_ = &newB;
  }

  double operator()(int i, int j) const {
    return (*operand_)(i,j) + u_[i];
  }
};

namespace pressio{
template<>
struct Traits<JacobianAction>{
  using scalar_type = double;
};

namespace ops{
void product(::pressio::transpose , ::pressio::nontranspose,
	     double alpha,
	     const Eigen::MatrixXd & A,
	     const JacobianAction & B,
	     double beta,
	     Eigen::MatrixXd & M)
{

  for (auto i=0; i<A.cols(); ++i){
    for (auto j=0; j<A.cols(); ++j){
      double sum=0.;
      for (auto k=0; k<A.rows(); ++k){
	sum += A(k,i) * B(k,j);
      }
      M(i,j) = beta * M(i,j) + alpha*sum;
    }
  }
}
}// end namespace ops
}// end namespace pressio

#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_steady.hpp"

struct MyFom
{
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

  state_type createState() const{ return state_type(N_); }
  residual_type createResidual() const{ return residual_type(N_); }

  JacobianAction createApplyJacobianResult(const Eigen::MatrixXd & B) const{
    return JacobianAction(N_, B);
  }

  void residual(const state_type & u, residual_type & r) const{
    EXPECT_TRUE(u.size()==r.size());
    EXPECT_TRUE(u.size()==N_);
    for (auto i=0; i<r.rows(); ++i){
     r(i) = u(i) + 1.;
    }
  }

  void applyJacobian(const state_type & state,
		     const Eigen::MatrixXd & B,
		     JacobianAction & A) const
  {
    A.store(state, B);
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
      system.residualAndJacobian(state, R, J);
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
      system.residualAndJacobian(state, R, J);
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

TEST(rom_galerkin_steady, test2)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  constexpr int N = 8;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  using full_state_type = typename fom_t::state_type;
  auto space = pressio::rom::create_trial_subspace<reduced_state_type, full_state_type>(phi);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = pressio::rom::galerkin::create_default_problem(space, fomSystem);

  FakeNonLinSolverSteady nonLinSolver(N);
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
