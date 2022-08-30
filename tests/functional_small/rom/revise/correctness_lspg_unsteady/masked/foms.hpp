
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_MASKED_CORRECT_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_MASKED_CORRECT_FOMS_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"

struct TrivialFomContTimeEigen
{
  using scalar_type    = double;
  using state_type     = Eigen::VectorXd;
  using velocity_type  = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomContTimeEigen(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
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
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomDiscreteTimeEigen(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     R(it) = -4354.323;
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
    // this is just faking things, to reproduce the cont-time above
    auto appJac = B;
    appJac.array() += time;
    A = (B - dt*appJac);
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(2434.3);
    }
  }
};


struct TrivialFomContTimeCustomTypes
{
  using scalar_type    = double;
  using state_type     = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type  = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomContTimeCustomTypes(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

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

    for (int i=0; i<N_; ++i){
     f(i) = u(i) + time;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
      for (std::size_t j=0; j< A.extent(1); ++j){
       A(it,j) = -1114;
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
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomDiscreteTimeCustomTypes(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

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

    for (std::size_t i=0; i<R.extent(0); ++i){
      auto f = y_np1(i) + time;
      R(i) = y_np1(i) -y_n(i) - dt*f;
    }

    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     R(it) = -4354.323;
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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
      for (std::size_t j=0; j< A.extent(1); ++j){
       A(it,j) = 2434.3;
     }
    }
  }
};

#endif
