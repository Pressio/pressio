
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_DEFAULT_CORRECT_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_DEFAULT_CORRECT_FOMS_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

struct TrivialFomContTimeEigen
{
  using scalar_type    = double;
  using state_type     = Eigen::VectorXd;
  using velocity_type  = state_type;
  int N_ = {};

  TrivialFomContTimeEigen(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)N_);

    for (int i=0; i<f.rows(); ++i){
     f(i) = u(i) + time;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     scalar_type time,
                     OperandType & A) const
  {
    A = B;
    A.array() += time;
  }
};


struct TrivialFomDiscreteTimeEigen
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using discrete_time_residual_type = state_type;
  using phi_type = Eigen::MatrixXd;
  int N_ = {};

  TrivialFomDiscreteTimeEigen(int N) : N_(N){}

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
    discrete_time_residual_type f(R.size());
    f.setZero();
    for (auto i=0; i<f.rows(); ++i){
     f(i) = y_np1(i) + time;
    }

    R = y_np1 -y_n - dt*f;
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
    // this is just faking things, to reproduce the cont-time above
    auto appJac = B;
    appJac.array() += time;
    A = (B - dt*appJac);
  }
};


struct TrivialFomContTimeCustomTypes
{
  using scalar_type    = double;
  using state_type     = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type  = state_type;
  int N_ = {};

  TrivialFomContTimeCustomTypes(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.extent(1));
    return A;
  }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.extent(0)==(std::size_t)f.extent(0));
    EXPECT_TRUE((std::size_t)u.extent(0)==(std::size_t)N_);

    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = u(i) + time;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     scalar_type time,
                     OperandType & A) const
  {
    A = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) += time;
      }
    }
  }
};


struct TrivialFomDiscreteTimeCustomTypes
{
  using scalar_type = double;
  using state_type     = ::pressiotests::MyCustomVector<scalar_type>;
  using discrete_time_residual_type = state_type;
  using phi_type = ::pressiotests::MyCustomMatrix<scalar_type>;
  int N_ = {};

  TrivialFomDiscreteTimeCustomTypes(int N) : N_(N){}

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
    discrete_time_residual_type f(R.extent(0));
    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = y_np1(i) + time;
    }

    for (std::size_t i=0; i<f.extent(0); ++i){
      R(i) = y_np1(i) -y_n(i) - dt*f(i);
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
    // this has no physical meaning, just to reproduce the cont-time above

    auto appJac = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        appJac(i,j) += time;
      }
    }

    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) = B(i,j) - dt*appJac(i,j);
      }
    }
  }
};

#endif
