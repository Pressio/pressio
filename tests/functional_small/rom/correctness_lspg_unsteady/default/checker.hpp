
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_DEFAULT_CORRECTNESS_CHECKER_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_DEFAULT_CORRECTNESS_CHECKER_HPP_

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


template<class phi_t>
struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};
  phi_t phi_;
  double dt_;

  FakeNonLinSolver(int N, phi_t phi, double dt)
    : N_(N), phi_(phi), dt_(dt){}

  template<class FomY_t, class RomStateType>
  void compute_y(FomY_t & y, const RomStateType & yrom)
  {
    for (int i=0; i<N_; ++i){
      y(i) = 0;
      for (int j=0; j<3; ++j){
	y(i) += phi_(i,j) * yrom(j);
      }
    }
  }

  template<class FomY_t, class F_t>
  void compute_f(const FomY_t & y, F_t & f, double time)
  {
    for (int i=0; i<N_; ++i){
      f(i) = y(i) + time;
    }
  }

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & romState)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);

    std::cout << " CALL COUNT = " << call_count_ << '\n';

    //*******************************
    //
    // call_count == 1
    //
    //*******************************
    if(call_count_==1)
    {
      double mytime = 2.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);
	// for callcount1, iter1: ynp1 and yn should be same, so use f only
	for (int i=0; i<N_; ++i){
	  EXPECT_DOUBLE_EQ(R(i), -dt_*f(i));
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_DOUBLE_EQ(J(i,j), phi_val - dt_*(phi_val+mytime));
	  }
	}

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
  // std::cout << "R = \n" << R << std::endl;
  // std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);

	StateType romState_n(3);
	romState_n(0) = 0; romState_n(1) = 1; romState_n(2) = 2;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, romState_n);

	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);

	for (int i=0; i<N_; ++i){
	  EXPECT_DOUBLE_EQ(R(i), ynp1(i)-yn(i)-dt_*f(i));
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_DOUBLE_EQ(J(i,j), phi_val - dt_*(phi_val+2.));
	  }
	}

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
      double mytime = 4.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
  // std::cout << "R = \n" << R << std::endl;
  // std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);
	// for callcount2, iter1: ynp1 and yn should be same, so use f only
	for (int i=0; i<N_; ++i){
	  EXPECT_DOUBLE_EQ(R(i), -dt_*f(i));
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_DOUBLE_EQ(J(i,j), phi_val - dt_*(phi_val+mytime));
	  }
	}

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residual(romState, R);
	system.jacobian(romState, J);
  // std::cout << "R = \n" << R << std::endl;
  // std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);

	StateType romState_n(3);
	romState_n(0) = 2; romState_n(1) = 3; romState_n(2) = 4.;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, romState_n);

	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);

	for (int i=0; i<N_; ++i){
	  EXPECT_DOUBLE_EQ(R(i), ynp1(i)-yn(i)-dt_*f(i));
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_DOUBLE_EQ(J(i,j), phi_val - dt_*(phi_val+mytime));
	  }
	}

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }
    }

  }
};

#endif
