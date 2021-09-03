
#ifndef PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_FOMS_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_FOMS_HPP_

#include <gtest/gtest.h>

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "pressio/ops.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#endif

struct TrivialFomOnlyVelocityCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomOnlyVelocityCustomTypes(int N, std::vector<int> ind)
  : N_(N), indices_to_corrupt_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const{
    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = u(i) + time;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
    }
  }
};

struct TrivialFomVelocityAndJacobianCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomVelocityAndJacobianCustomTypes(int N, std::vector<int> ind)
  : N_(N), indices_to_corrupt_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class scalar_type>
  pressiotests::MyCustomMatrix<scalar_type> 
  createApplyJacobianResult(const pressiotests::MyCustomMatrix<scalar_type> & B) const
  {
    pressiotests::MyCustomMatrix<scalar_type> A(B.extent(0), B.extent(1));
    return A;
  }

  // computes: A = Jac B
  void applyJacobian(const state_type & state,
                     const pressiotests::MyCustomMatrix<scalar_type> & B,
                     const scalar_type & time,
                     pressiotests::MyCustomMatrix<scalar_type> & A) const
  {
    A = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) += time;
      }
    }
    for (auto & it : indices_to_corrupt_){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(it,j) = -4232;
      }
    }
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = u(i) + time;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
    }
  }
};


struct TrivialFomDiscreteTimeCustomTypes
{
  using scalar_type = double;
  using state_type = ::pressiotests::MyCustomVector<scalar_type>;
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
    EXPECT_TRUE((std::size_t)y_np1.extent(0)==(std::size_t)y_n.extent(0));
    EXPECT_TRUE((std::size_t)y_np1.extent(0)==(std::size_t)R.extent(0));
    EXPECT_TRUE((std::size_t)R.extent(0)==(std::size_t)N_);

    for (int i=0; i<N_; ++i){
     auto f = y_np1(i) + time;
     R(i) = y_np1(i) -y_n(i) - dt*f;
    }

    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     R(it) = -1114;
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
    EXPECT_TRUE((std::size_t)y_np1.extent(0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)y_n.extent(0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.extent(0)==(std::size_t)N_);

    A = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) += time;
      }
    }
    for (auto & it : indices_to_corrupt_){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(it,j) = -4232;
      }
    }
  }
};


struct TrivialFomOnlyVelocityEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomOnlyVelocityEigen(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const{
    for (auto i=0; i<f.rows(); ++i){
     f(i) = u(i) + time;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
    }
  }
};

struct TrivialFomVelocityAndJacobianEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomVelocityAndJacobianEigen(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    return OperandType(B.rows(), B.cols());
  }

  // computes: A = Jac B
  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     const scalar_type & time,
                     OperandType & A) const
  {
    A = B;
    A.array() += time;
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
    }
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (auto i=0; i<f.rows(); ++i){
     f(i) = u(i) + time;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
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
    EXPECT_TRUE((std::size_t)y_np1.size()==(std::size_t)y_n.size());
    EXPECT_TRUE((std::size_t)y_np1.size()==(std::size_t)R.size());
    EXPECT_TRUE((std::size_t)R.size()==(std::size_t)N_);

    for (int i=0; i<N_; ++i)
    {
     auto f = y_np1(i) + time;
     R(i) = y_np1(i) -y_n(i) - dt*f;
    }

    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     R(it) = -1114;
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
    EXPECT_TRUE((std::size_t)y_np1.size()==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)y_n.size()==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    A = B;
    A.array() += time;
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
    }
  }
};


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
struct TrivialFomOnlyVelocityTpetra
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
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomOnlyVelocityTpetra(int N, std::vector<int> ind)
  : N_(N), indices_to_corrupt_(ind)
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

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    f.putScalar(time);
    pressio::ops::update(f, 1, u, 1);

    auto lv = f.getLocalViewHost();
    auto mygids = contigMap_->getMyGlobalIndices();
    for (std::size_t i=0; i<mygids.extent(0); ++i){
      if (std::find(indices_to_corrupt_.cbegin(), indices_to_corrupt_.cend(), mygids(i))!=indices_to_corrupt_.cend()){
        lv(i,0) = 241.;
      }
    }
  }
};
#endif

#endif
