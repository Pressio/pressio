
#ifndef PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_COMMON_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_COMMON_HPP_

#include <gtest/gtest.h>

#include "../custom_data_types.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "pressio/ops.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif


// ==============================================
// ==============================================
//
//  FOM CLASSES
//
// ==============================================
// ==============================================

struct TrivialFomOnlyVelocityEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomOnlyVelocityEigen(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==indices.size());
  }

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + time;
    }
  }
};

struct TrivialFomVelocityAndJacobianEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomVelocityAndJacobianEigen(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  // computes: A = Jac B
  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     const scalar_type & time,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)state.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
    A.array() += time;
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t)N_);
    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + time;
    }
  }
};

struct TrivialFomDiscreteTimeEigen
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using discrete_time_residual_type = state_type;
  using phi_type = Eigen::MatrixXd;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomDiscreteTimeEigen(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  discrete_time_residual_type createDiscreteTimeResidual() const{ 
    return discrete_time_residual_type(N_);
  }

  phi_type createApplyDiscreteTimeJacobianResult(const phi_type & B) const{ 
    return phi_type(N_, B.cols());
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type & R, 
                              const state_type & y_np1,
                              const state_type & y_n) const
  {
    EXPECT_TRUE((std::size_t)y_np1.size()==(std::size_t)y_n.size());
    EXPECT_TRUE((std::size_t)y_np1.size()!=(std::size_t)R.size());
    EXPECT_TRUE((std::size_t)R.size()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i)
    {
     auto f = y_np1(indices[i]) + time;
     R(i) = y_np1(indices[i]) -y_n(indices[i]) - dt*f;
    }
  }

  template<class StepCountType>
  void applyDiscreteTimeJacobian(StepCountType, 
                              double time, 
                              double dt, 
                              const phi_type & B, 
                              phi_type & A,
                              const state_type & y_np1,
                              const state_type & y_n) const
  {
    EXPECT_TRUE((std::size_t)y_np1.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)y_n.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
    A.array() += time;
  }
};


struct TrivialFomOnlyVelocityCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomOnlyVelocityCustomTypes(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + time;
    }
  }
};

struct TrivialFomVelocityAndJacobianCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomVelocityAndJacobianCustomTypes(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class scalar_type>
  pressiotests::MyCustomMatrix<scalar_type> 
  createApplyJacobianResult(const pressiotests::MyCustomMatrix<scalar_type> & B) const
  {
    pressiotests::MyCustomMatrix<scalar_type> A(N_, B.extent(1));
    return A;
  }

  // computes: A = Jac B
  void applyJacobian(const state_type & state,
                     const pressiotests::MyCustomMatrix<scalar_type> & B,
                     const scalar_type & time,
                     pressiotests::MyCustomMatrix<scalar_type> & A) const
  {
    for (std::size_t i=0; i<indices.size(); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) = B(indices[i], j) + time;
      }
    }
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + time;
    }
  }
};

struct TrivialFomDiscreteTimeCustomTypes
{
  using scalar_type = double;
  using state_type = ::pressiotests::MyCustomVector<scalar_type>;
  using discrete_time_residual_type = state_type;
  using phi_type = pressiotests::MyCustomMatrix<scalar_type>;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomDiscreteTimeCustomTypes(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  discrete_time_residual_type createDiscreteTimeResidual() const{ 
    return discrete_time_residual_type(N_);
  }

  phi_type createApplyDiscreteTimeJacobianResult(const phi_type & B) const{ 
    return phi_type(N_, B.extent(1));
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type & R, 
                              const state_type & y_np1,
                              const state_type & y_n) const
  {
    EXPECT_TRUE((std::size_t)y_np1.extent(0)==(std::size_t)y_n.extent(0));
    EXPECT_TRUE((std::size_t)y_np1.extent(0)!=(std::size_t)R.extent(0));
    EXPECT_TRUE((std::size_t)R.extent(0)==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i)
    {
     auto f = y_np1(indices[i]) + time;
     R(i) = y_np1(indices[i]) -y_n(indices[i]) - dt*f;
    }
  }

  template<class StepCountType>
  void applyDiscreteTimeJacobian(StepCountType, 
                              double time, 
                              double dt, 
                              const phi_type & B, 
                              phi_type & A,
                              const state_type & y_np1,
                              const state_type & y_n) const
  {
    EXPECT_TRUE((std::size_t)y_np1.extent(0)!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)y_n.extent(0)!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.extent(0)==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) = B(indices[i], j) + time;
      }
    }
  }
};

// ==============================================
// ==============================================
//
//  PROJECTORS
//
// ==============================================
// ==============================================

// for explicit, projector acts on FOM velicity
struct ProjectorExplicitEigen
{
  using operator_type = Eigen::MatrixXd; 
  operator_type matrix_;

  ProjectorExplicitEigen(const operator_type & phi) : matrix_(phi){}

  template<class operand_type, class ScalarType, class ResultType>
  void operator()(const operand_type & operand, ScalarType time, ResultType & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

using ProjectorImplicitEigen = ProjectorExplicitEigen;


// for explicit, projector acts on FOM velicity
template<class ScalarType>
struct ProjectorExplicitCustomTypes
{
  using operator_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  operator_type matrix_;

  ProjectorExplicitCustomTypes(const operator_type & phi) : matrix_(phi){}

  // result is the projected RHS, so it is a rom type
  template<class operand_type>
  void operator()(const operand_type & operand, ScalarType time, Eigen::VectorXd & result) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t k=0; k<matrix_.extent(1); ++k)
    {
      result(k) = 0;
      for (std::size_t i=0; i<matrix_.extent(0); ++i){
        result(k) += matrix_(i,k)*operand(i);
      }
    }
  }
};

// for explicit, projector acts on FOM velicity
template<class ScalarType>
struct ProjectorImplicitCustomTypes
{
  using operator_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  operator_type matrix_;

  ProjectorImplicitCustomTypes(const operator_type & phi) : matrix_(phi){}

  void operator()(const ::pressiotests::MyCustomVector<ScalarType> & operand, 
             ScalarType time, 
             Eigen::VectorXd & result) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t k=0; k<matrix_.extent(1); ++k)
    {
      result(k) = 0;
      for (std::size_t i=0; i<matrix_.extent(0); ++i){
        result(k) += matrix_(i,k)*operand(i);
      }
    }
  }

  void operator()(const ::pressiotests::MyCustomMatrix<ScalarType> & operand, 
             ScalarType time, 
             Eigen::MatrixXd & result) const
  {
    for (std::size_t i=0; i<matrix_.extent(1); ++i){
      for (std::size_t j=0; j<operand.extent(1); ++j)
      {
        result(i,j) = 0;
        for (std::size_t k=0; k<matrix_.extent(0); ++k){
          result(i,j) += matrix_(k,i)*operand(k,j);
        }
      }
    }
  }
};

struct ObserverA
{
  void operator()(int32_t step, double time, Eigen::VectorXd state)
  {
    EXPECT_TRUE(step<=2);

    if (step==0){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 1.);
      EXPECT_DOUBLE_EQ(state[2], 2.);
    }
    if (step==1){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 51.);
      EXPECT_DOUBLE_EQ(state[2], 102.);
    }
    if (step==2){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 2611.);
      EXPECT_DOUBLE_EQ(state[2], 5222.);
    }
  }
};

struct FakeNonLinSolverContTime
{
  int call_count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << R << std::endl;
      // std::cout << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 0.);
      EXPECT_DOUBLE_EQ(R[1], -140.);
      EXPECT_DOUBLE_EQ(R[2], -280.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << R << std::endl;
      // std::cout << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 1.);
      EXPECT_DOUBLE_EQ(R[1], -199.);
      EXPECT_DOUBLE_EQ(R[2], -399.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }

    //
    // call_count == 2
    //
    if(call_count_==2)
    {
      // do solver iterator 1
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << R << std::endl;
      // std::cout << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0],    0.);
      EXPECT_DOUBLE_EQ(R[1], -300.);
      EXPECT_DOUBLE_EQ(R[2], -600.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << R << std::endl;
      // std::cout << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0],    1.);
      EXPECT_DOUBLE_EQ(R[1], -359.);
      EXPECT_DOUBLE_EQ(R[2], -719.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

struct FakeNonLinSolverForDiscreteTime
{
  int call_count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 0.);
      EXPECT_DOUBLE_EQ(R[1], -140.);
      EXPECT_DOUBLE_EQ(R[2], -280.);

      EXPECT_DOUBLE_EQ(J(0,0),  0.);
      EXPECT_DOUBLE_EQ(J(0,1),  0.);
      EXPECT_DOUBLE_EQ(J(0,2),  0.);
      EXPECT_DOUBLE_EQ(J(1,0), 20.);
      EXPECT_DOUBLE_EQ(J(1,1), 30.);
      EXPECT_DOUBLE_EQ(J(1,2), 40.);
      EXPECT_DOUBLE_EQ(J(2,0), 40.);
      EXPECT_DOUBLE_EQ(J(2,1), 60.);
      EXPECT_DOUBLE_EQ(J(2,2), 80.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 0.);
      EXPECT_DOUBLE_EQ(R[1], -170.);
      EXPECT_DOUBLE_EQ(R[2], -340.);

      EXPECT_DOUBLE_EQ(J(0,0),  0.);
      EXPECT_DOUBLE_EQ(J(0,1),  0.);
      EXPECT_DOUBLE_EQ(J(0,2),  0.);
      EXPECT_DOUBLE_EQ(J(1,0), 20.);
      EXPECT_DOUBLE_EQ(J(1,1), 30.);
      EXPECT_DOUBLE_EQ(J(1,2), 40.);
      EXPECT_DOUBLE_EQ(J(2,0), 40.);
      EXPECT_DOUBLE_EQ(J(2,1), 60.);
      EXPECT_DOUBLE_EQ(J(2,2), 80.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }

    //
    // call_count == 2
    //
    if(call_count_==2)
    {
      // do solver iterator 1
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0],    0.);
      EXPECT_DOUBLE_EQ(R[1], -300.);
      EXPECT_DOUBLE_EQ(R[2], -600.);

      EXPECT_DOUBLE_EQ(J(0,0),  0.);
      EXPECT_DOUBLE_EQ(J(0,1),  0.);
      EXPECT_DOUBLE_EQ(J(0,2),  0.);
      EXPECT_DOUBLE_EQ(J(1,0), 40.);
      EXPECT_DOUBLE_EQ(J(1,1), 50.);
      EXPECT_DOUBLE_EQ(J(1,2), 60.);
      EXPECT_DOUBLE_EQ(J(2,0), 80.);
      EXPECT_DOUBLE_EQ(J(2,1), 100.);
      EXPECT_DOUBLE_EQ(J(2,2), 120.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residual(state, R);
      system.jacobian(state, J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0],    0.);
      EXPECT_DOUBLE_EQ(R[1], -330.);
      EXPECT_DOUBLE_EQ(R[2], -660.);

      EXPECT_DOUBLE_EQ(J(0,0),  0.);
      EXPECT_DOUBLE_EQ(J(0,1),  0.);
      EXPECT_DOUBLE_EQ(J(0,2),  0.);
      EXPECT_DOUBLE_EQ(J(1,0), 40.);
      EXPECT_DOUBLE_EQ(J(1,1), 50.);
      EXPECT_DOUBLE_EQ(J(1,2), 60.);
      EXPECT_DOUBLE_EQ(J(2,0), 80.);
      EXPECT_DOUBLE_EQ(J(2,1), 100.);
      EXPECT_DOUBLE_EQ(J(2,2), 120.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};

#endif
