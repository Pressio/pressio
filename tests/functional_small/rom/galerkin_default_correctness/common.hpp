
#ifndef PRESSIO_TEST_ROM_GALERKIN_DEFAULT_CORRECT_COMMON_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_DEFAULT_CORRECT_COMMON_HPP_

#include <gtest/gtest.h>

#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

struct TrivialFomOnlyVelocityCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};

  TrivialFomOnlyVelocityCustomTypes(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const{
    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = u(i) + time;
    }
  }
};

struct TrivialFomVelocityAndJacobianCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using velocity_type     = state_type;
  int N_ = {};

  TrivialFomVelocityAndJacobianCustomTypes(int N): N_(N){}

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
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (std::size_t i=0; i<f.extent(0); ++i){
     f(i) = u(i) + time;
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
     R(i) = y_np1(i) - y_n(i) - dt*f(i);
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
    A = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) += time;
      }
    }
  }
};


struct TrivialFomOnlyVelocityEigen
{
  using scalar_type	      = double;
  using state_type	      = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};

  TrivialFomOnlyVelocityEigen(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const{
    for (auto i=0; i<f.rows(); ++i){
	   f(i) = u(i) + time;
    }
  }
};

struct TrivialFomVelocityAndJacobianEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using velocity_type     = state_type;
  int N_ = {};

  TrivialFomVelocityAndJacobianEigen(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(B.rows(), B.cols());
    return A;
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
  }

  void velocity(const state_type & u, const scalar_type time, velocity_type & f) const
  {
    for (auto i=0; i<f.rows(); ++i){
     f(i) = u(i) + time;
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
    A = B;
    A.array() += time;
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

  TrivialFomOnlyVelocityTpetra(int N): N_(N)
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
  }
};

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
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
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
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
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
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
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
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
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
