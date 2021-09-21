
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_MASKED_CORRECT_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_MASKED_CORRECT_FOMS_HPP_

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
  const std::vector<int> indices_ = {};

  TrivialFomContTimeEigen(std::vector<int> ind)
    : N_(ind.size()), indices_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t)N_);
    for (std::size_t i=0; i<indices_.size(); ++i){
     f(i) = u(indices_[i]) + time;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     scalar_type time,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)state.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    for (std::size_t i=0; i<indices_.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices_[i], j) + time;
      }
    }
  }
};

struct TrivialFomContTimeCustomTypes
{
  using scalar_type    = double;
  using state_type     = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type  = state_type;
  int N_ = {};
  const std::vector<int> indices_ = {};

  TrivialFomContTimeCustomTypes(std::vector<int> ind)
    : N_(ind.size()), indices_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.extent(1));
    return A;
  }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<indices_.size(); ++i){
     f(i) = u(indices_[i]) + time;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     scalar_type time,
                     OperandType & A) const
  {
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) = B(indices_[i],j) + time;
      }
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
  const std::vector<int> indices_ = {};

  TrivialFomDiscreteTimeEigen(std::vector<int> ind)
    : N_(ind.size()), indices_(ind){}

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
    for (std::size_t i=0; i<indices_.size(); ++i)
    {
      auto f = y_np1(indices_[i]) + time;
      R(i) = y_np1(indices_[i]) -y_n(indices_[i]) - dt*f;
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
    // compute J*B where J is the FOM jacobian
    // see cont-time adapter above for what J*B is
    phi_type JB(N_, B.cols());
    for (std::size_t i=0; i< (std::size_t)N_; ++i){
      for (std::size_t j=0; j< (std::size_t)JB.cols(); ++j){
        JB(i,j) = B(indices_[i],j) + time;
      }
    }

    // apply discrete-time Jacobian to B
    // thisis like: A = I - dt*JB since we are mimicing BDF1
    for (std::size_t i=0; i<indices_.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices_[i], j) - dt *JB(i,j);
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
  const std::vector<int> indices_ = {};

  TrivialFomDiscreteTimeCustomTypes(std::vector<int> ind)
    : N_(ind.size()), indices_(ind){}

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

    for (std::size_t i=0; i<indices_.size(); ++i)
    {
      auto f = y_np1(indices_[i]) + time;
      R(i) = y_np1(indices_[i]) -y_n(indices_[i]) - dt*f;
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
    // compute J*B where J is the FOM jacobian
    // see cont-time adapter above for what J*B is
    phi_type JB(N_, B.extent(1));
    for (int i=0; i< N_; ++i){
      for (std::size_t j=0; j< JB.extent(1); ++j){
        JB(i,j) = B(indices_[i],j) + time;
      }
    }

    // apply discrete-time Jacobian to B
    // thisis like: A = I - dt*JB since we are mimicing BDF1
    for (std::size_t i=0; i<indices_.size(); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) = B(indices_[i], j) - dt *JB(i,j);
      }
    }
  }
};

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
struct TrivialFomVelocityAndJacobianTpetra
{
  using scalar_type       = double;
  using state_type        = Tpetra::Vector<>;
  using velocity_type     = state_type;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;
  int N_ = {};
  int rank_;
  int numProc_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;

  TrivialFomVelocityAndJacobianTpetra(int N): N_(N)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    contigMap_ = Teuchos::rcp(new map_t(N_, 0, comm_));
  }

  Teuchos::RCP<const tcomm> comm(){ return comm_; }
  Teuchos::RCP<const map_t> map(){ return contigMap_; }

  velocity_type createVelocity() const
  {
    vec_t result(contigMap_);
    return result;
  }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(contigMap_, B.getNumVectors());
    return A;
  }

  // computes: A = Jac B
  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     const scalar_type & time,
                     OperandType & A) const
  {
    auto tmp = pressio::ops::clone(B);
    tmp.putScalar(time);

    pressio::ops::deep_copy(A,B);
    pressio::ops::update(A,1,tmp,1);
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    f.putScalar(time);
    pressio::ops::update(f, 1, u, 1);
  }
};
#endif

#endif
