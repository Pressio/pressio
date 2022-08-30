
#ifndef PRESSIO_TEST_ROM_GALERKIN_CORRECTNESS_CHECKERS_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_CORRECTNESS_CHECKERS_HPP_

#include <gtest/gtest.h>
#include "pressio/ode.hpp"

// struct FakeNonLinSolverContTime
// {
//   int call_count_ = 0;

//   template<class SystemType, class StateType>
//   void solve(const SystemType & system, StateType & state)
//   {
//     ++call_count_;
//     auto R = system.createResidual();
//     auto J = system.createJacobian();

//     //
//     // call_count == 1
//     //
//     if(call_count_==1)
//     {
//       // do solver iterator 1
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0], 0.);
//       EXPECT_DOUBLE_EQ(R[1], -140.);
//       EXPECT_DOUBLE_EQ(R[2], -280.);

//       EXPECT_DOUBLE_EQ(J(0,0),   1.);
//       EXPECT_DOUBLE_EQ(J(0,1),   0.);
//       EXPECT_DOUBLE_EQ(J(0,2),   0.);
//       EXPECT_DOUBLE_EQ(J(1,0), -40.);
//       EXPECT_DOUBLE_EQ(J(1,1), -59.);
//       EXPECT_DOUBLE_EQ(J(1,2), -80.);
//       EXPECT_DOUBLE_EQ(J(2,0), -80.);
//       EXPECT_DOUBLE_EQ(J(2,1),-120.);
//       EXPECT_DOUBLE_EQ(J(2,2),-159.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

//       // do solver iterator 2
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0], 1.);
//       EXPECT_DOUBLE_EQ(R[1], -199.);
//       EXPECT_DOUBLE_EQ(R[2], -399.);

//       EXPECT_DOUBLE_EQ(J(0,0),   1.);
//       EXPECT_DOUBLE_EQ(J(0,1),   0.);
//       EXPECT_DOUBLE_EQ(J(0,2),   0.);
//       EXPECT_DOUBLE_EQ(J(1,0), -40.);
//       EXPECT_DOUBLE_EQ(J(1,1), -59.);
//       EXPECT_DOUBLE_EQ(J(1,2), -80.);
//       EXPECT_DOUBLE_EQ(J(2,0), -80.);
//       EXPECT_DOUBLE_EQ(J(2,1),-120.);
//       EXPECT_DOUBLE_EQ(J(2,2),-159.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
//     }

//     //
//     // call_count == 2
//     //
//     if(call_count_==2)
//     {
//       // do solver iterator 1
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0],    0.);
//       EXPECT_DOUBLE_EQ(R[1], -300.);
//       EXPECT_DOUBLE_EQ(R[2], -600.);

//       EXPECT_DOUBLE_EQ(J(0,0),   1.);
//       EXPECT_DOUBLE_EQ(J(0,1),   0.);
//       EXPECT_DOUBLE_EQ(J(0,2),   0.);
//       EXPECT_DOUBLE_EQ(J(1,0), -80.);
//       EXPECT_DOUBLE_EQ(J(1,1), -99.);
//       EXPECT_DOUBLE_EQ(J(1,2), -120.);
//       EXPECT_DOUBLE_EQ(J(2,0), -160.);
//       EXPECT_DOUBLE_EQ(J(2,1), -200.);
//       EXPECT_DOUBLE_EQ(J(2,2), -239.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

//       // do solver iterator 2
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0],    1.);
//       EXPECT_DOUBLE_EQ(R[1], -359.);
//       EXPECT_DOUBLE_EQ(R[2], -719.);

//       EXPECT_DOUBLE_EQ(J(0,0),   1.);
//       EXPECT_DOUBLE_EQ(J(0,1),   0.);
//       EXPECT_DOUBLE_EQ(J(0,2),   0.);
//       EXPECT_DOUBLE_EQ(J(1,0), -80.);
//       EXPECT_DOUBLE_EQ(J(1,1), -99.);
//       EXPECT_DOUBLE_EQ(J(1,2), -120.);
//       EXPECT_DOUBLE_EQ(J(2,0), -160.);
//       EXPECT_DOUBLE_EQ(J(2,1), -200.);
//       EXPECT_DOUBLE_EQ(J(2,2), -239.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
//     }
//   }
// };


// struct FakeNonLinSolverForDiscreteTime
// {
//   int call_count_ = 0;

//   template<class SystemType, class StateType>
//   void solve(const SystemType & system, StateType & state)
//   {
//     ++call_count_;
//     auto R = system.createResidual();
//     auto J = system.createJacobian();

//     //
//     // call_count == 1
//     //
//     if(call_count_==1)
//     {
//       // do solver iterator 1
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0], 0.);
//       EXPECT_DOUBLE_EQ(R[1], -140.);
//       EXPECT_DOUBLE_EQ(R[2], -280.);

//       EXPECT_DOUBLE_EQ(J(0,0),  0.);
//       EXPECT_DOUBLE_EQ(J(0,1),  0.);
//       EXPECT_DOUBLE_EQ(J(0,2),  0.);
//       EXPECT_DOUBLE_EQ(J(1,0), 20.);
//       EXPECT_DOUBLE_EQ(J(1,1), 30.);
//       EXPECT_DOUBLE_EQ(J(1,2), 40.);
//       EXPECT_DOUBLE_EQ(J(2,0), 40.);
//       EXPECT_DOUBLE_EQ(J(2,1), 60.);
//       EXPECT_DOUBLE_EQ(J(2,2), 80.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

//       // do solver iterator 2
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0], 0.);
//       EXPECT_DOUBLE_EQ(R[1], -170.);
//       EXPECT_DOUBLE_EQ(R[2], -340.);

//       EXPECT_DOUBLE_EQ(J(0,0),  0.);
//       EXPECT_DOUBLE_EQ(J(0,1),  0.);
//       EXPECT_DOUBLE_EQ(J(0,2),  0.);
//       EXPECT_DOUBLE_EQ(J(1,0), 20.);
//       EXPECT_DOUBLE_EQ(J(1,1), 30.);
//       EXPECT_DOUBLE_EQ(J(1,2), 40.);
//       EXPECT_DOUBLE_EQ(J(2,0), 40.);
//       EXPECT_DOUBLE_EQ(J(2,1), 60.);
//       EXPECT_DOUBLE_EQ(J(2,2), 80.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
//     }

//     //
//     // call_count == 2
//     //
//     if(call_count_==2)
//     {
//       // do solver iterator 1
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0],    0.);
//       EXPECT_DOUBLE_EQ(R[1], -300.);
//       EXPECT_DOUBLE_EQ(R[2], -600.);

//       EXPECT_DOUBLE_EQ(J(0,0),  0.);
//       EXPECT_DOUBLE_EQ(J(0,1),  0.);
//       EXPECT_DOUBLE_EQ(J(0,2),  0.);
//       EXPECT_DOUBLE_EQ(J(1,0), 40.);
//       EXPECT_DOUBLE_EQ(J(1,1), 50.);
//       EXPECT_DOUBLE_EQ(J(1,2), 60.);
//       EXPECT_DOUBLE_EQ(J(2,0), 80.);
//       EXPECT_DOUBLE_EQ(J(2,1), 100.);
//       EXPECT_DOUBLE_EQ(J(2,2), 120.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

//       // do solver iterator 2
//       system.residual(state, R);
//       system.jacobian(state, J);
//       // std::cout << "S " << call_count_ << " \n" << R << std::endl;
//       // std::cout << "S " << call_count_ << " \n" << J << std::endl;
//       EXPECT_DOUBLE_EQ(R[0],    0.);
//       EXPECT_DOUBLE_EQ(R[1], -330.);
//       EXPECT_DOUBLE_EQ(R[2], -660.);

//       EXPECT_DOUBLE_EQ(J(0,0),  0.);
//       EXPECT_DOUBLE_EQ(J(0,1),  0.);
//       EXPECT_DOUBLE_EQ(J(0,2),  0.);
//       EXPECT_DOUBLE_EQ(J(1,0), 40.);
//       EXPECT_DOUBLE_EQ(J(1,1), 50.);
//       EXPECT_DOUBLE_EQ(J(1,2), 60.);
//       EXPECT_DOUBLE_EQ(J(2,0), 80.);
//       EXPECT_DOUBLE_EQ(J(2,1), 100.);
//       EXPECT_DOUBLE_EQ(J(2,2), 120.);

//       for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
//     }
//   }
// };

#endif
