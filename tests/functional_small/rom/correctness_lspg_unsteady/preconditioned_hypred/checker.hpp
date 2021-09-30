
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_HYPRED_CORRECTNESS_CHECKER_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_HYPRED_CORRECTNESS_CHECKER_HPP_

#include <gtest/gtest.h>

struct ObserverA
{
  template<class StateType>
  void operator()(int32_t step, double time, const StateType & state)
  {
    EXPECT_TRUE(step<=2);

    if (step==0){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 1.);
      EXPECT_DOUBLE_EQ(state[2], 2.);
    }
    if (step==1){
      EXPECT_DOUBLE_EQ(state[0], 2.);
      EXPECT_DOUBLE_EQ(state[1], 3.);
      EXPECT_DOUBLE_EQ(state[2], 4.);
    }
    if (step==2){
      EXPECT_DOUBLE_EQ(state[0], 4.);
      EXPECT_DOUBLE_EQ(state[1], 5.);
      EXPECT_DOUBLE_EQ(state[2], 6.);
    }
  }
};

struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};

  FakeNonLinSolver(int N) : N_(N){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & romState)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    //*******************************
    //
    // call_count == 1
    //
    //*******************************
    if(call_count_==1)
    {
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -14.+1.5);
	EXPECT_DOUBLE_EQ(R(1), -32.+1.5);
	EXPECT_DOUBLE_EQ(R(2), -50.+1.5);
	EXPECT_DOUBLE_EQ(R(3), -68.+1.5);
	EXPECT_DOUBLE_EQ(R(4), -86.+1.5);
	EXPECT_DOUBLE_EQ(R(5), -104.+1.5);
	EXPECT_DOUBLE_EQ(R(6), -122.+1.5);
	EXPECT_DOUBLE_EQ(R(7), -140.+1.5);

	EXPECT_DOUBLE_EQ(J(0,0), -4.+1.5);
	EXPECT_DOUBLE_EQ(J(0,1), -5.+1.5);
	EXPECT_DOUBLE_EQ(J(0,2), -6.+1.5);
	EXPECT_DOUBLE_EQ(J(1,0), -7.+1.5);
	EXPECT_DOUBLE_EQ(J(1,1), -8.+1.5);
	EXPECT_DOUBLE_EQ(J(1,2), -9.+1.5);
	EXPECT_DOUBLE_EQ(J(2,0), -10.+1.5);
	EXPECT_DOUBLE_EQ(J(2,1), -11.+1.5);
	EXPECT_DOUBLE_EQ(J(2,2), -12.+1.5);
	EXPECT_DOUBLE_EQ(J(3,0), -13.+1.5);
	EXPECT_DOUBLE_EQ(J(3,1), -14.+1.5);
	EXPECT_DOUBLE_EQ(J(3,2), -15.+1.5);
	EXPECT_DOUBLE_EQ(J(4,0), -16.+1.5);
	EXPECT_DOUBLE_EQ(J(4,1), -17.+1.5);
	EXPECT_DOUBLE_EQ(J(4,2), -18.+1.5);
	EXPECT_DOUBLE_EQ(J(5,0), -19.+1.5);
	EXPECT_DOUBLE_EQ(J(5,1), -20.+1.5);
	EXPECT_DOUBLE_EQ(J(5,2), -21.+1.5);
	EXPECT_DOUBLE_EQ(J(6,0), -22.+1.5);
	EXPECT_DOUBLE_EQ(J(6,1), -23.+1.5);
	EXPECT_DOUBLE_EQ(J(6,2), -24.+1.5);
	EXPECT_DOUBLE_EQ(J(7,0), -25.+1.5);
	EXPECT_DOUBLE_EQ(J(7,1), -26.+1.5);
	EXPECT_DOUBLE_EQ(J(7,2), -27.+1.5);

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -17.+1.5);
	EXPECT_DOUBLE_EQ(R(1), -44.+1.5);
	EXPECT_DOUBLE_EQ(R(2), -71.+1.5);
	EXPECT_DOUBLE_EQ(R(3), -98.+1.5);
	EXPECT_DOUBLE_EQ(R(4), -125.+1.5);
	EXPECT_DOUBLE_EQ(R(5), -152.+1.5);
	EXPECT_DOUBLE_EQ(R(6), -179.+1.5);
	EXPECT_DOUBLE_EQ(R(7), -206.+1.5);

	EXPECT_DOUBLE_EQ(J(0,0), -4.+1.5);
	EXPECT_DOUBLE_EQ(J(0,1), -5.+1.5);
	EXPECT_DOUBLE_EQ(J(0,2), -6.+1.5);
	EXPECT_DOUBLE_EQ(J(1,0), -7.+1.5);
	EXPECT_DOUBLE_EQ(J(1,1), -8.+1.5);
	EXPECT_DOUBLE_EQ(J(1,2), -9.+1.5);
	EXPECT_DOUBLE_EQ(J(2,0), -10.+1.5);
	EXPECT_DOUBLE_EQ(J(2,1), -11.+1.5);
	EXPECT_DOUBLE_EQ(J(2,2), -12.+1.5);
	EXPECT_DOUBLE_EQ(J(3,0), -13.+1.5);
	EXPECT_DOUBLE_EQ(J(3,1), -14.+1.5);
	EXPECT_DOUBLE_EQ(J(3,2), -15.+1.5);
	EXPECT_DOUBLE_EQ(J(4,0), -16.+1.5);
	EXPECT_DOUBLE_EQ(J(4,1), -17.+1.5);
	EXPECT_DOUBLE_EQ(J(4,2), -18.+1.5);
	EXPECT_DOUBLE_EQ(J(5,0), -19.+1.5);
	EXPECT_DOUBLE_EQ(J(5,1), -20.+1.5);
	EXPECT_DOUBLE_EQ(J(5,2), -21.+1.5);
	EXPECT_DOUBLE_EQ(J(6,0), -22.+1.5);
	EXPECT_DOUBLE_EQ(J(6,1), -23.+1.5);
	EXPECT_DOUBLE_EQ(J(6,2), -24.+1.5);
	EXPECT_DOUBLE_EQ(J(7,0), -25.+1.5);
	EXPECT_DOUBLE_EQ(J(7,1), -26.+1.5);
	EXPECT_DOUBLE_EQ(J(7,2), -27.+1.5);

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }
    }

    //*******************************
    //
    // call_count == 2
    //
    //*******************************
    if(call_count_==2)
    {
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -30.+1.5);
	EXPECT_DOUBLE_EQ(R(1), -84.+1.5);
	EXPECT_DOUBLE_EQ(R(2), -138.+1.5);
	EXPECT_DOUBLE_EQ(R(3), -192.+1.5);
	EXPECT_DOUBLE_EQ(R(4), -246.+1.5);
	EXPECT_DOUBLE_EQ(R(5), -300.+1.5);
	EXPECT_DOUBLE_EQ(R(6), -354.+1.5);
	EXPECT_DOUBLE_EQ(R(7), -408.+1.5);

	EXPECT_DOUBLE_EQ(J(0,0), -8.+1.5);
	EXPECT_DOUBLE_EQ(J(0,1), -9.+1.5);
	EXPECT_DOUBLE_EQ(J(0,2), -10.+1.5);
	EXPECT_DOUBLE_EQ(J(1,0), -11.+1.5);
	EXPECT_DOUBLE_EQ(J(1,1), -12.+1.5);
	EXPECT_DOUBLE_EQ(J(1,2), -13.+1.5);
	EXPECT_DOUBLE_EQ(J(2,0), -14.+1.5);
	EXPECT_DOUBLE_EQ(J(2,1), -15.+1.5);
	EXPECT_DOUBLE_EQ(J(2,2), -16.+1.5);
	EXPECT_DOUBLE_EQ(J(3,0), -17.+1.5);
	EXPECT_DOUBLE_EQ(J(3,1), -18.+1.5);
	EXPECT_DOUBLE_EQ(J(3,2), -19.+1.5);
	EXPECT_DOUBLE_EQ(J(4,0), -20.+1.5);
	EXPECT_DOUBLE_EQ(J(4,1), -21.+1.5);
	EXPECT_DOUBLE_EQ(J(4,2), -22.+1.5);
	EXPECT_DOUBLE_EQ(J(5,0), -23.+1.5);
	EXPECT_DOUBLE_EQ(J(5,1), -24.+1.5);
	EXPECT_DOUBLE_EQ(J(5,2), -25.+1.5);
	EXPECT_DOUBLE_EQ(J(6,0), -26.+1.5);
	EXPECT_DOUBLE_EQ(J(6,1), -27.+1.5);
	EXPECT_DOUBLE_EQ(J(6,2), -28.+1.5);
	EXPECT_DOUBLE_EQ(J(7,0), -29.+1.5);
	EXPECT_DOUBLE_EQ(J(7,1), -30.+1.5);
	EXPECT_DOUBLE_EQ(J(7,2), -31.+1.5);

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -33.+1.5);
	EXPECT_DOUBLE_EQ(R(1), -96.+1.5);
	EXPECT_DOUBLE_EQ(R(2), -159.+1.5);
	EXPECT_DOUBLE_EQ(R(3), -222.+1.5);
	EXPECT_DOUBLE_EQ(R(4), -285.+1.5);
	EXPECT_DOUBLE_EQ(R(5), -348.+1.5);
	EXPECT_DOUBLE_EQ(R(6), -411.+1.5);
	EXPECT_DOUBLE_EQ(R(7), -474.+1.5);

	EXPECT_DOUBLE_EQ(J(0,0), -8.+1.5);
	EXPECT_DOUBLE_EQ(J(0,1), -9.+1.5);
	EXPECT_DOUBLE_EQ(J(0,2), -10.+1.5);
	EXPECT_DOUBLE_EQ(J(1,0), -11.+1.5);
	EXPECT_DOUBLE_EQ(J(1,1), -12.+1.5);
	EXPECT_DOUBLE_EQ(J(1,2), -13.+1.5);
	EXPECT_DOUBLE_EQ(J(2,0), -14.+1.5);
	EXPECT_DOUBLE_EQ(J(2,1), -15.+1.5);
	EXPECT_DOUBLE_EQ(J(2,2), -16.+1.5);
	EXPECT_DOUBLE_EQ(J(3,0), -17.+1.5);
	EXPECT_DOUBLE_EQ(J(3,1), -18.+1.5);
	EXPECT_DOUBLE_EQ(J(3,2), -19.+1.5);
	EXPECT_DOUBLE_EQ(J(4,0), -20.+1.5);
	EXPECT_DOUBLE_EQ(J(4,1), -21.+1.5);
	EXPECT_DOUBLE_EQ(J(4,2), -22.+1.5);
	EXPECT_DOUBLE_EQ(J(5,0), -23.+1.5);
	EXPECT_DOUBLE_EQ(J(5,1), -24.+1.5);
	EXPECT_DOUBLE_EQ(J(5,2), -25.+1.5);
	EXPECT_DOUBLE_EQ(J(6,0), -26.+1.5);
	EXPECT_DOUBLE_EQ(J(6,1), -27.+1.5);
	EXPECT_DOUBLE_EQ(J(6,2), -28.+1.5);
	EXPECT_DOUBLE_EQ(J(7,0), -29.+1.5);
	EXPECT_DOUBLE_EQ(J(7,1), -30.+1.5);
	EXPECT_DOUBLE_EQ(J(7,2), -31.+1.5);

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }
    }

  }
};


struct FakeNonLinSolverTpetra
{
  int call_count_ = 0;

  FakeNonLinSolverTpetra(){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & romState)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    // EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    // EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    // EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    //*******************************
    //
    // call_count == 1
    //
    //*******************************
    if(call_count_==1)
    {
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;
	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }
    }

    //*******************************
    //
    // call_count == 2
    //
    //*******************************
    if(call_count_==2)
    {
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;
	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;
	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }
    }

  }
};

#endif
