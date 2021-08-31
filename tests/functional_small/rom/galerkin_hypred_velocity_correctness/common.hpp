
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

struct TrivialFomOnlyVelocityEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomOnlyVelocityEigen(int N) : N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
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

  TrivialFomVelocityAndJacobianEigen(int N) : N_(N){}

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
    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
    A.array() += time;
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + time;
    }
  }
};

// for explicit, projector acts on FOM velicity
struct ProjectorExplicitEigen
{
  using operator_type = Eigen::MatrixXd; 
  operator_type matrix_;

  ProjectorExplicitEigen(const operator_type & phi) : matrix_(phi){}

  template<class operand_type, class ResultType>
  void apply(const operand_type & operand, ResultType & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

using ProjectorImplicitEigen = ProjectorExplicitEigen;


struct TrivialFomOnlyVelocityCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  TrivialFomOnlyVelocityCustomTypes(int N) : N_(N){}

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

  TrivialFomVelocityAndJacobianCustomTypes(int N) : N_(N){}

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

// for explicit, projector acts on FOM velicity
template<class ScalarType>
struct ProjectorExplicitCustomTypes
{
  using operator_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  operator_type matrix_;

  ProjectorExplicitCustomTypes(const operator_type & phi) : matrix_(phi){}

  // result is the projected RHS, so it is a rom type
  template<class operand_type>
  void apply(const operand_type & operand, Eigen::VectorXd & result) const
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

  void apply(const ::pressiotests::MyCustomVector<ScalarType> & operand, 
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

  void apply(const ::pressiotests::MyCustomMatrix<ScalarType> & operand, 
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

struct FakeNonLinSolver
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

#endif
