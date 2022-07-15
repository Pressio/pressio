
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

	EXPECT_DOUBLE_EQ(R(0), -14.);
	EXPECT_DOUBLE_EQ(R(1), -32.);
	EXPECT_DOUBLE_EQ(R(2), -50.);
	EXPECT_DOUBLE_EQ(R(3), -68.);
	EXPECT_DOUBLE_EQ(R(4), -86.);
	EXPECT_DOUBLE_EQ(R(5), -104.);
	EXPECT_DOUBLE_EQ(R(6), -122.);
	EXPECT_DOUBLE_EQ(R(7), -140.);

	EXPECT_DOUBLE_EQ(J(0,0), -4.);
	EXPECT_DOUBLE_EQ(J(0,1), -5.);
	EXPECT_DOUBLE_EQ(J(0,2), -6.);
	EXPECT_DOUBLE_EQ(J(1,0), -7.);
	EXPECT_DOUBLE_EQ(J(1,1), -8.);
	EXPECT_DOUBLE_EQ(J(1,2), -9.);
	EXPECT_DOUBLE_EQ(J(2,0), -10.);
	EXPECT_DOUBLE_EQ(J(2,1), -11.);
	EXPECT_DOUBLE_EQ(J(2,2), -12.);
	EXPECT_DOUBLE_EQ(J(3,0), -13.);
	EXPECT_DOUBLE_EQ(J(3,1), -14.);
	EXPECT_DOUBLE_EQ(J(3,2), -15.);
	EXPECT_DOUBLE_EQ(J(4,0), -16.);
	EXPECT_DOUBLE_EQ(J(4,1), -17.);
	EXPECT_DOUBLE_EQ(J(4,2), -18.);
	EXPECT_DOUBLE_EQ(J(5,0), -19.);
	EXPECT_DOUBLE_EQ(J(5,1), -20.);
	EXPECT_DOUBLE_EQ(J(5,2), -21.);
	EXPECT_DOUBLE_EQ(J(6,0), -22.);
	EXPECT_DOUBLE_EQ(J(6,1), -23.);
	EXPECT_DOUBLE_EQ(J(6,2), -24.);
	EXPECT_DOUBLE_EQ(J(7,0), -25.);
	EXPECT_DOUBLE_EQ(J(7,1), -26.);
	EXPECT_DOUBLE_EQ(J(7,2), -27.);

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

	EXPECT_DOUBLE_EQ(R(0), -17.);
	EXPECT_DOUBLE_EQ(R(1), -44.);
	EXPECT_DOUBLE_EQ(R(2), -71.);
	EXPECT_DOUBLE_EQ(R(3), -98.);
	EXPECT_DOUBLE_EQ(R(4), -125.);
	EXPECT_DOUBLE_EQ(R(5), -152.);
	EXPECT_DOUBLE_EQ(R(6), -179.);
	EXPECT_DOUBLE_EQ(R(7), -206.);

	EXPECT_DOUBLE_EQ(J(0,0), -4.);
	EXPECT_DOUBLE_EQ(J(0,1), -5.);
	EXPECT_DOUBLE_EQ(J(0,2), -6.);
	EXPECT_DOUBLE_EQ(J(1,0), -7.);
	EXPECT_DOUBLE_EQ(J(1,1), -8.);
	EXPECT_DOUBLE_EQ(J(1,2), -9.);
	EXPECT_DOUBLE_EQ(J(2,0), -10.);
	EXPECT_DOUBLE_EQ(J(2,1), -11.);
	EXPECT_DOUBLE_EQ(J(2,2), -12.);
	EXPECT_DOUBLE_EQ(J(3,0), -13.);
	EXPECT_DOUBLE_EQ(J(3,1), -14.);
	EXPECT_DOUBLE_EQ(J(3,2), -15.);
	EXPECT_DOUBLE_EQ(J(4,0), -16.);
	EXPECT_DOUBLE_EQ(J(4,1), -17.);
	EXPECT_DOUBLE_EQ(J(4,2), -18.);
	EXPECT_DOUBLE_EQ(J(5,0), -19.);
	EXPECT_DOUBLE_EQ(J(5,1), -20.);
	EXPECT_DOUBLE_EQ(J(5,2), -21.);
	EXPECT_DOUBLE_EQ(J(6,0), -22.);
	EXPECT_DOUBLE_EQ(J(6,1), -23.);
	EXPECT_DOUBLE_EQ(J(6,2), -24.);
	EXPECT_DOUBLE_EQ(J(7,0), -25.);
	EXPECT_DOUBLE_EQ(J(7,1), -26.);
	EXPECT_DOUBLE_EQ(J(7,2), -27.);

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

	EXPECT_DOUBLE_EQ(R(0), -30.);
	EXPECT_DOUBLE_EQ(R(1), -84.);
	EXPECT_DOUBLE_EQ(R(2), -138.);
	EXPECT_DOUBLE_EQ(R(3), -192.);
	EXPECT_DOUBLE_EQ(R(4), -246.);
	EXPECT_DOUBLE_EQ(R(5), -300.);
	EXPECT_DOUBLE_EQ(R(6), -354.);
	EXPECT_DOUBLE_EQ(R(7), -408.);

	EXPECT_DOUBLE_EQ(J(0,0), -8.);
	EXPECT_DOUBLE_EQ(J(0,1), -9.);
	EXPECT_DOUBLE_EQ(J(0,2), -10.);
	EXPECT_DOUBLE_EQ(J(1,0), -11.);
	EXPECT_DOUBLE_EQ(J(1,1), -12.);
	EXPECT_DOUBLE_EQ(J(1,2), -13.);
	EXPECT_DOUBLE_EQ(J(2,0), -14.);
	EXPECT_DOUBLE_EQ(J(2,1), -15.);
	EXPECT_DOUBLE_EQ(J(2,2), -16.);
	EXPECT_DOUBLE_EQ(J(3,0), -17.);
	EXPECT_DOUBLE_EQ(J(3,1), -18.);
	EXPECT_DOUBLE_EQ(J(3,2), -19.);
	EXPECT_DOUBLE_EQ(J(4,0), -20.);
	EXPECT_DOUBLE_EQ(J(4,1), -21.);
	EXPECT_DOUBLE_EQ(J(4,2), -22.);
	EXPECT_DOUBLE_EQ(J(5,0), -23.);
	EXPECT_DOUBLE_EQ(J(5,1), -24.);
	EXPECT_DOUBLE_EQ(J(5,2), -25.);
	EXPECT_DOUBLE_EQ(J(6,0), -26.);
	EXPECT_DOUBLE_EQ(J(6,1), -27.);
	EXPECT_DOUBLE_EQ(J(6,2), -28.);
	EXPECT_DOUBLE_EQ(J(7,0), -29.);
	EXPECT_DOUBLE_EQ(J(7,1), -30.);
	EXPECT_DOUBLE_EQ(J(7,2), -31.);

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

	EXPECT_DOUBLE_EQ(R(0), -33.);
	EXPECT_DOUBLE_EQ(R(1), -96.);
	EXPECT_DOUBLE_EQ(R(2), -159.);
	EXPECT_DOUBLE_EQ(R(3), -222.);
	EXPECT_DOUBLE_EQ(R(4), -285.);
	EXPECT_DOUBLE_EQ(R(5), -348.);
	EXPECT_DOUBLE_EQ(R(6), -411.);
	EXPECT_DOUBLE_EQ(R(7), -474.);

	EXPECT_DOUBLE_EQ(J(0,0), -8.);
	EXPECT_DOUBLE_EQ(J(0,1), -9.);
	EXPECT_DOUBLE_EQ(J(0,2), -10.);
	EXPECT_DOUBLE_EQ(J(1,0), -11.);
	EXPECT_DOUBLE_EQ(J(1,1), -12.);
	EXPECT_DOUBLE_EQ(J(1,2), -13.);
	EXPECT_DOUBLE_EQ(J(2,0), -14.);
	EXPECT_DOUBLE_EQ(J(2,1), -15.);
	EXPECT_DOUBLE_EQ(J(2,2), -16.);
	EXPECT_DOUBLE_EQ(J(3,0), -17.);
	EXPECT_DOUBLE_EQ(J(3,1), -18.);
	EXPECT_DOUBLE_EQ(J(3,2), -19.);
	EXPECT_DOUBLE_EQ(J(4,0), -20.);
	EXPECT_DOUBLE_EQ(J(4,1), -21.);
	EXPECT_DOUBLE_EQ(J(4,2), -22.);
	EXPECT_DOUBLE_EQ(J(5,0), -23.);
	EXPECT_DOUBLE_EQ(J(5,1), -24.);
	EXPECT_DOUBLE_EQ(J(5,2), -25.);
	EXPECT_DOUBLE_EQ(J(6,0), -26.);
	EXPECT_DOUBLE_EQ(J(6,1), -27.);
	EXPECT_DOUBLE_EQ(J(6,2), -28.);
	EXPECT_DOUBLE_EQ(J(7,0), -29.);
	EXPECT_DOUBLE_EQ(J(7,1), -30.);
	EXPECT_DOUBLE_EQ(J(7,2), -31.);

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
