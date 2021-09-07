
#ifndef PRESSIO_TEST_ROM_LSPG_STEADY_CORRECTNESS_CHECKER_HPP_
#define PRESSIO_TEST_ROM_LSPG_STEADY_CORRECTNESS_CHECKER_HPP_

#include <gtest/gtest.h>

struct ObserverA
{
  void operator()(int32_t step, double time, Eigen::VectorXd state)
  {
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

struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};

  FakeNonLinSolver(int N) : N_(N){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE(pressio::ops::extent(R,0)==N_);
    EXPECT_TRUE(pressio::ops::extent(J,0)==N_);
    EXPECT_TRUE(pressio::ops::extent(J,1)==3);
    std::cout << "S " << call_count_ << " \n" << R << std::endl;

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residual(state, R);
      system.jacobian(state, J);
      std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      // EXPECT_DOUBLE_EQ(R(0), 0*0.+1*1.+2*2.+1.);
      // EXPECT_DOUBLE_EQ(R(1), 3*0.+4*1.+5*2.+1.);
      // EXPECT_DOUBLE_EQ(R(2), 6*0.+7*1.+8*2.+1.);
      // EXPECT_DOUBLE_EQ(R(3), 9*0.+10*1.+11*2.+1.);
      // EXPECT_DOUBLE_EQ(R(4), 12*0.+13*1.+14*2.+1.);
      // EXPECT_DOUBLE_EQ(R(5), 15*0.+16*1.+17*2.+1.);
      // EXPECT_DOUBLE_EQ(R(6), 18*0.+19*1.+20*2.+1.);
      // EXPECT_DOUBLE_EQ(R(7), 21*0.+22*1.+23*2.+1.);

      // double start = 1;
      // int count = 0;
      // for (int i=0; i<N_; ++i){
      //   for (int j=0; j<3; ++j){
      //     EXPECT_DOUBLE_EQ(J(i,j), start + (double)count++);
      //   }
      // }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2 
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      // EXPECT_DOUBLE_EQ(R(0), 0*1.+1*2.+2*3.+1.);
      // EXPECT_DOUBLE_EQ(R(1), 3*1.+4*2.+5*3.+1.);
      // EXPECT_DOUBLE_EQ(R(2), 6*1.+7*2.+8*3.+1.);
      // EXPECT_DOUBLE_EQ(R(3), 9*1.+10*2.+11*3.+1.);
      // EXPECT_DOUBLE_EQ(R(4), 12*1.+13*2.+14*3.+1.);
      // EXPECT_DOUBLE_EQ(R(5), 15*1.+16*2.+17*3.+1.);
      // EXPECT_DOUBLE_EQ(R(6), 18*1.+19*2.+20*3.+1.);
      // EXPECT_DOUBLE_EQ(R(7), 21*1.+22*2.+23*3.+1.);

      // start = 1;
      // count = 0;
      // for (int i=0; i<N_; ++i){
      //   for (int j=0; j<3; ++j){
      //     EXPECT_DOUBLE_EQ(J(i,j), start + (double)count++);
      //   }
      // }

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

#endif
