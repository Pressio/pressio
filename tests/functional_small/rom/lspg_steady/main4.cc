
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_steady.hpp"

struct MyFom
{
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  MyFom(int N,  std::vector<int> ind): N_(N), indices_to_corrupt_(ind){}

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

  void applyJacobian(const state_type & /*unused*/,
                     const Eigen::MatrixXd & B,
                     Eigen::MatrixXd & A) const{
    A = B;
    for (int i=0; i<A.rows(); ++i){
      for (int j=0; j<A.cols(); ++j){
        A(i,j) += 1.;
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
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, J, true);

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
      system.residualAndJacobian(state, R, J, true);

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

    }// end if call_count == 1
  }// end solve
};

class MyMasker
{
  const std::vector<int> sample_indices_ = {};

public:
  MyMasker(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  auto createApplyMaskResult(const Eigen::VectorXd & /*operand*/) const{
    return Eigen::VectorXd(sample_indices_.size());
  }

  auto createApplyMaskResult(const Eigen::MatrixXd & operand) const{
    return Eigen::MatrixXd(sample_indices_.size(), operand.cols());
  }

  template<class T1, class T2>
  void operator()(const Eigen::MatrixBase<T1> & operand,
		  Eigen::MatrixBase<T2> & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.cols(); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};


// struct MyMaskerResidual
// {
//   const std::vector<int> sample_indices_ = {};
//   using operand_type = Eigen::VectorXd;
//   using result_type = operand_type;

//   MyMaskerResidual(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

//   result_type createApplyMaskResult(const operand_type & /*operand*/) const{
//     return result_type(sample_indices_.size());
//   }

//   void operator()(const operand_type & operand, result_type & result) const
//   {
//     for (std::size_t i=0; i<sample_indices_.size(); ++i){
//       result(i) = operand(sample_indices_[i]);
//     }
//   }
// };

// struct MyMaskerJacAction
// {
//   const std::vector<int> sample_indices_ = {};
//   using operand_type = Eigen::MatrixXd;
//   using result_type = operand_type;

//   MyMaskerJacAction(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

//   result_type createApplyMaskResult(const operand_type & operand) const{
//     return result_type(sample_indices_.size(), operand.cols());
//   }

//   void operator()(const operand_type & operand, result_type & result) const
//   {
//     for (std::size_t i=0; i<sample_indices_.size(); ++i){
//       for (int j=0; j<operand.cols(); ++j){
//         result(i,j) = operand(sample_indices_[i],j);
//       }
//     }
//   }
// };

TEST(rom_lspg_steady, test4)
{
  /* steady masked LSPG */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  const int N = 15;
  const std::vector<int> indices_to_corrupt = {1,3,5,7,9,11,13};

  using fom_t = MyFom;
  fom_t fomSystem(N, indices_to_corrupt);

  Eigen::MatrixXd phi(N, 3);
  int count = 0;
  for (int i=0; i<phi.rows(); ++i){
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
  typename fom_t::state_type shift(N);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};
  MyMasker  masker(sample_indices);
  auto problem = pressio::rom::lspg::create_steady_problem(space, fomSystem, masker);

  FakeNonLinSolverSteady nonLinSolver(sample_indices.size());
  nonLinSolver.solve(problem, romState);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 2.);
  EXPECT_DOUBLE_EQ(romState[1], 3.);
  EXPECT_DOUBLE_EQ(romState[2], 4.);

  pressio::log::finalize();
}
