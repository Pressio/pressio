
#ifndef PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_FOMS_COMMON_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_FOMS_COMMON_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "pressio/ops.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

struct TrivialFomOnlyVelocityEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};

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

#endif
