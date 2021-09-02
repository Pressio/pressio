
#ifndef PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_HPP_

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
#include <Tpetra_Import_decl.hpp>
#endif

// ==============================================
// ==============================================
//
//  FOM CLASSES
//
// ==============================================
// ==============================================

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
    EXPECT_TRUE(y_np1.extent(0)==y_n.extent(0));
    EXPECT_TRUE(y_np1.extent(0)==R.extent(0));
    EXPECT_TRUE(R.extent(0)==N_);

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
    EXPECT_TRUE(y_np1.extent(0)==N_);
    EXPECT_TRUE(y_n.extent(0)==N_);
    EXPECT_TRUE(A.extent(0)==N_);

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
    EXPECT_TRUE(y_np1.size()==y_n.size());
    EXPECT_TRUE(y_np1.size()==R.size());
    EXPECT_TRUE(R.size()==N_);

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
    EXPECT_TRUE(y_np1.size()==N_);
    EXPECT_TRUE(y_n.size()==N_);
    EXPECT_TRUE(A.rows()==N_);

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



// ==============================================
// ==============================================
//
//  MASKERS 
//
// ==============================================
// ==============================================

// for explicit, masker acts on FOM velicity
struct MaskerExplicitEigen
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = typename TrivialFomOnlyVelocityEigen::velocity_type;

  MaskerExplicitEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(sample_indices_.size());
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

// for implicit, masker acts on FOM velicity and FOM apply jac result
struct MaskerImplicitEigen
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  MaskerImplicitEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.cols());
  }

  template<class TimeType>
  void operator()(const vec_operand_type & operand, TimeType time, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  template<class TimeType>
  void operator()(const mat_operand_type & operand, TimeType time, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.cols(); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

// for explicit, masker acts on FOM velicity
template<class ScalarType>
struct MaskerExplicitCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = ::pressiotests::MyCustomVector<ScalarType>;

  MaskerExplicitCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(sample_indices_.size());
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

// for implicit, masker acts on FOM velicity and FOM apply jac result
template<class ScalarType>
struct MaskerImplicitCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = ::pressiotests::MyCustomVector<ScalarType>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<ScalarType>;

  MaskerImplicitCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.extent(1));
  }

  template<class TimeType>
  void operator()(const vec_operand_type & operand, TimeType time, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  template<class TimeType>
  void operator()(const mat_operand_type & operand, TimeType time, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (std::size_t j=0; j<operand.extent(1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// for explicit, masker acts on FOM velicity
struct MaskerExplicitTpetra
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using GO = typename map_t::global_ordinal_type;
  using LO = typename map_t::local_ordinal_type;

  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> maskedMap_;
  std::vector<int> sample_indices_ = {};
  using operand_type = typename TrivialFomOnlyVelocityTpetra::velocity_type;

  MaskerExplicitTpetra(Teuchos::RCP<const tcomm> comm, 
                       Teuchos::RCP<const map_t> full_map,                       
                       std::vector<int> sample_indices) 
  : comm_(comm), sample_indices_(sample_indices)
  {
    auto mygids_fm = full_map->getMyGlobalIndices();
    // 100 is large enough for this test, only a subset 
    // will be picked by the map constructor
    GO mymaskedGids[100];
    GO count=0;
    for (std::size_t i=0; i<mygids_fm.extent(0); ++i)
    {
      const bool r = std::find(sample_indices_.cbegin(), sample_indices_.cend(), mygids_fm(i))!=sample_indices_.cend();
      if (r){ mymaskedGids[count++] = mygids_fm(i); }
    }

    maskedMap_ = Teuchos::rcp(new map_t(-1, mymaskedGids, count, 0, comm_));
  }

  Teuchos::RCP<const map_t> map(){ return maskedMap_; }

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(maskedMap_);
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    using importer_t = Tpetra::Import<LO, GO>;
    importer_t importer(operand.getMap(), maskedMap_);
    result.doImport(operand, importer, Tpetra::CombineMode::REPLACE);
  }
};
#endif




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

  template<class operand_type, class ScalarType,  class ResultType>
  void operator()(const operand_type & operand, const ScalarType time, ResultType & result) const
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// for explicit, projector acts on FOM velicity
struct ProjectorExplicitTpetra
{
  using operator_type = Tpetra::MultiVector<>; 
  operator_type matrix_;

  ProjectorExplicitTpetra(const operator_type & phi) : matrix_(phi){}

  template<class operand_type, class ScalarType, class ResultType>
  void operator()(const operand_type & operand, const ScalarType time, ResultType & result) const
  {
    pressio::ops::product(pressio::transpose{}, 1, matrix_, operand, 0, result);
  }
};
#endif



// ==============================================
// ==============================================
//
//  REST
//
// ==============================================
// ==============================================

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
